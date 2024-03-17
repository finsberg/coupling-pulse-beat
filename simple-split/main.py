"""Attempt to solve mechanics and ep models separately and then couple them together using the missing values
similar to the approach in the 0D model. However, since we use different meshes for the mechanics and electrophysiology
we will need to interpolate the missing values between the meshes."""

import gotranx
from typing import Any
from dataclasses import dataclass
import numpy as np
import dolfin
import beat
import pulse
import logging
import ufl_legacy as ufl
import matplotlib.pyplot as plt


import mechanicssolver


@dataclass
class MissingValue:
    element: ufl.FiniteElementBase
    mechanics_mesh: dolfin.Mesh
    ep_mesh: dolfin.Mesh

    def __post_init__(self):
        self.V_ep = dolfin.FunctionSpace(self.ep_mesh, self.element)
        self.V_mechanics = dolfin.FunctionSpace(self.mechanics_mesh, self.element)

        self.u_ep = dolfin.Function(self.V_ep)
        self.u_mechanics = dolfin.Function(self.V_mechanics)

        self.values_ep = np.zeros((1, self.u_ep.vector().local_size()))
        self.values_mechanics = np.zeros((1, self.u_mechanics.vector().local_size()))

    def ep_values_to_function(self) -> None:
        self.u_ep.vector().set_local(self.values_ep[0])

    def ep_function_to_values(self) -> None:
        self.values_ep[0, :] = self.u_ep.vector().get_local()

    def mechanics_values_to_function(self) -> None:
        self.u_mechanics.vector().set_local(self.values_mechanics[0])

    def mechanics_function_to_values(self) -> None:
        self.values_mechanics[0, :] = self.u_mechanics.vector().get_local()

    def interpolate_ep_to_mechanics(self) -> None:
        self.u_mechanics.interpolate(self.u_ep)

    def interpolate_mechanics_to_ep(self) -> None:
        self.u_ep.interpolate(self.u_mechanics)

    @property
    def ep_array(self):
        return self.u_ep.vector().get_local()

    @ep_array.setter
    def ep_array(self, value):
        self.u_ep.vector().set_local(value)

    @property
    def mechanics_array(self):
        return self.u_mechanics.vector().get_local()

    @mechanics_array.setter
    def mechanics_array(self, value):
        self.u_mechanics.vector().set_local(value)


def setup_geometry(dx):
    Lx = 20.0  # mm
    Ly = 7.0  # mm
    Lz = 3.0  # mm

    mesh = dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
        int(np.rint((Lz / dx))),
    )
    return mesh


def define_conductivity_tensor(chi, C_m):
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = 0.17  # mS / mm
    sigma_it = 0.019  # mS / mm
    sigma_el = 0.62  # mS / mm
    sigma_et = 0.24  # mS / mm

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(C_m * chi)
    s_l = sigma_l / (C_m * chi)  # mm^2 / ms
    s_t = sigma_t / (C_m * chi)  # mm^2 / ms

    # Define conductivity tensor
    M = dolfin.as_tensor(((s_l, 0, 0), (0, s_t, 0), (0, 0, s_t)))

    return M


def define_stimulus(mesh, chi, C_m, time, A=50000.0, duration=2.0):
    S1_marker = 1
    L = 1.5
    S1_subdomain = dolfin.CompiledSubDomain(
        "x[0] <= L + DOLFIN_EPS && x[1] <= L + DOLFIN_EPS && x[2] <= L + DOLFIN_EPS",
        L=L,
    )
    S1_markers = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    S1_subdomain.mark(S1_markers, S1_marker)

    # Define stimulation (NB: region of interest carried by the mesh
    # and assumptions in cbcbeat)
    # duration = 2.0  # ms
    # A = 50000.0  # mu A/cm^3
    cm2mm = 10.0
    factor = 1.0 / (chi * C_m)  # NB: cbcbeat convention
    amplitude = factor * A * (1.0 / cm2mm) ** 3  # mV/ms

    I_s = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=0.0,
        duration=duration,
        amplitude=amplitude,
        degree=0,
    )

    dx = dolfin.Measure("dx", domain=mesh, subdomain_data=S1_markers)(S1_marker)
    return beat.base_model.Stimulus(dz=dx, expr=I_s)


# Load the model
ode = gotranx.load_ode("ORdmm_Land.ode")

mechanics_comp = ode.get_component("mechanics")
mechanics_ode = mechanics_comp.to_ode()

ep_ode = ode - mechanics_comp


# Generate code for the electrophysiology model
code_ep = gotranx.cli.gotran2py.get_code(
    ep_ode,
    scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen],
    missing_values=mechanics_ode.missing_variables,
)
# Generate code for the mechanics model
code_mechanics = gotranx.cli.gotran2py.get_code(
    mechanics_ode,
    scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen],
    missing_values=ep_ode.missing_variables,
)

# # Execute code and get the models
ep_model: dict[str, Any] = {}
exec(code_ep, ep_model)
mechanics_model: dict[str, Any] = {}
exec(code_mechanics, mechanics_model)


# Set time step to 0.1 ms
dt = 0.1
# Simulate model for 1000 ms
t = np.arange(0, 1000, dt)

# Get the index of the membrane potential
V_index_ep = ep_model["state_index"]("v")
# Forwared generalized rush larsen scheme for the electrophysiology model
fgr_ep = ep_model["forward_generalized_rush_larsen"]
# Monitor function for the electrophysiology model
mon_ep = ep_model["monitor"]
# Missing values function for the electrophysiology model
mv_ep = ep_model["missing_values"]
# Index of the calcium concentration
Ca_index_ep = ep_model["state_index"]("cai")

# Forwared generalized rush larsen scheme for the mechanics model
fgr_mechanics = mechanics_model["forward_generalized_rush_larsen"]
# Monitor function for the mechanics model
mon_mechanics = mechanics_model["monitor"]
# Missing values function for the mechanics model
mv_mechanics = mechanics_model["missing_values"]
# Index of the active tension
Ta_index_mechanics = mechanics_model["monitor_index"]("Ta")
# Index of the J_TRPN
J_TRPN_index_mechanics = mechanics_model["monitor_index"]("J_TRPN")


tol = 2e-5
# Create arrays to store the results
V_ep = np.zeros(len(t))
Ca_ep = np.zeros(len(t))

Ta_mechanics = np.zeros(len(t))
J_TRPN_mechanics = np.zeros(len(t))

fig, ax = plt.subplots(3, 2, sharex=True, figsize=(10, 10))

# Get initial values from the EP model
y_ep_ = ep_model["init_state_values"]()
p_ep = ep_model["init_parameter_values"](amp=0.0)
ep_missing_values_ = np.zeros(len(ep_ode.missing_variables))

# Get initial values from the mechanics model
y_mechanics_ = mechanics_model["init_state_values"]()
p_mechanics = mechanics_model["init_parameter_values"]()
mechanics_missing_values_ = np.zeros(len(mechanics_ode.missing_variables))


mesh = setup_geometry(dx=3.0)
ep_mesh = dolfin.adapt(dolfin.adapt(dolfin.adapt(mesh)))

# Surface to volume ratio
chi = 140.0  # mm^{-1}
# Membrane capacitance
C_m = 0.01  # mu F / mm^2

time = dolfin.Constant(0.0)
I_s = define_stimulus(mesh=ep_mesh, chi=chi, C_m=C_m, time=time, A=0)

M = define_conductivity_tensor(chi, C_m)

# params = {"linear_solver_type": "direct"}
params = {"preconditioner": "sor", "use_custom_preconditioner": False}
ep_ode_space = dolfin.FunctionSpace(ep_mesh, "Lagrange", 1)
v_ode = dolfin.Function(ep_ode_space)
num_points_ep = v_ode.vector().local_size()

mechanics_ode_space = dolfin.FunctionSpace(mesh, "Lagrange", 1)
v_ode_mechanics = dolfin.Function(mechanics_ode_space)
num_points_mechanics = v_ode_mechanics.vector().local_size()


missing_mech = MissingValue(
    element=mechanics_ode_space.ufl_element(),
    mechanics_mesh=mesh,
    ep_mesh=ep_mesh,
)

missing_ep = MissingValue(
    element=ep_ode_space.ufl_element(),
    mechanics_mesh=mesh,
    ep_mesh=ep_mesh,
)

# if len(mechanics_ode.missing_variables) == 1:
#     missing_mechanics_mech = dolfin.Function(mechanics_ode_space)
#     missing_mechanics_ep = dolfin.Function(ep_ode_space)
# else:
#     raise NotImplementedError("Only one missing variable is supported")

# if len(ep_ode.missing_variables) == 1:
#     missing_ep_ep = dolfin.Function(ep_ode_space)
#     missing_ep_mech = dolfin.Function(mechanics_ode_space)
# else:
#     raise NotImplementedError("Only one missing variable is supported")


y_ep = np.zeros((len(y_ep_), num_points_ep))
y_ep.T[:] = y_ep_

y_mechanics = np.zeros((len(y_mechanics_), num_points_mechanics))
y_mechanics.T[:] = y_mechanics_


# Get the default values of the missing values
# A little bit chicken and egg problem here, but in this specific case we know that
# the mechanics_missing_values is only the calcium concentration, which is a state variable
# and this doesn't require any additional information to be calculated.
mechanics_missing_values_[:] = mv_ep(0, y_ep_, p_ep, ep_missing_values_)
ep_missing_values_[:] = mv_mechanics(
    0, y_mechanics_, p_mechanics, mechanics_missing_values_
)

# ep_ep_missing_values = np.zeros((len(ep_ode.missing_variables), num_points_ep))
# ep_mech_missing_values = np.zeros((len(ep_ode.missing_variables), num_points_mechanics))
# mechanics_mech_missing_values = np.zeros(
#     (len(mechanics_ode.missing_variables), num_points_mechanics)
# )
# mechanics_ep_missing_values = np.zeros(
#     (len(mechanics_ode.missing_variables), num_points_ep)
# )

missing_mech.values_ep.T[:] = mechanics_missing_values_
missing_ep.values_mechanics.T[:] = ep_missing_values_
missing_ep.values_ep.T[:] = ep_missing_values_
missing_mech.values_mechanics.T[:] = mechanics_missing_values_
# mechanics_mech_missing_values.T[:] = mechanics_missing_values_
# mechanics_ep_missing_values.T[:] = mechanics_missing_values_

# ep_ep_missing_values.T[:] = ep_missing_values_
# ep_mech_missing_values.T[:] = ep_missing_values_

# We will store the previous missing values to check for convergence
# prev_mechanics_missing_values = np.zeros_like(mechanics_mech_missing_values)
# prev_mechanics_missing_values[:] = mechanics_mech_missing_values


pde = beat.MonodomainModel(time=time, mesh=ep_mesh, M=M, I_s=I_s, params=params)
ode = beat.odesolver.DolfinODESolver(
    v_ode=dolfin.Function(ep_ode_space),
    v_pde=pde.state,
    fun=fgr_ep,
    init_states=y_ep,
    parameters=p_ep,
    num_states=len(y_ep),
    v_index=ep_model["state_index"]("v"),
    missing_variables=missing_ep.values_ep,
    num_missing_variables=len(ep_ode.missing_variables),
)

ep_solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)


fixed = dolfin.CompiledSubDomain("on_boundary && near(x[0], 0.0)")
ffun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
ffun.set_all(0)
fixed.mark(ffun, 1)
marker_functions = pulse.MarkerFunctions(ffun=ffun)
f0 = dolfin.as_vector([1.0, 0.0, 0.0])
s0 = dolfin.as_vector([0.0, 1.0, 0.0])
n0 = dolfin.as_vector([0.0, 0.0, 1.0])
microstructure = pulse.Microstructure(f0=f0, s0=s0, n0=n0)
geometry = pulse.Geometry(
    mesh=mesh, marker_functions=marker_functions, microstructure=microstructure
)
# Use the default material parameters
material_parameters = pulse.HolzapfelOgden.default_parameters()

# Select model for active contraction
# active_model = pulse.ActiveModels.active_strain
active_model = "active_stress"

# Set the activation
V = dolfin.FunctionSpace(geometry.mesh, "Lagrange", 1)
activation = dolfin.Function(V)

# Create material
material = pulse.HolzapfelOgden(
    active_model=active_model,
    parameters=material_parameters,
    activation=activation,
)


# Make Dirichlet boundary conditions
def dirichlet_bc(W):
    V = W if W.sub(0).num_sub_spaces() == 0 else W.sub(0)
    return dolfin.DirichletBC(V, dolfin.Constant((0.0, 0.0, 0.0)), fixed)


# Collect Boundary Conditions
bcs = pulse.BoundaryConditions(dirichlet=(dirichlet_bc,))


# Create problem
mechanics_odesolver = mechanicssolver.ODESolver(
    Ta=activation,
    fun=fgr_mechanics,
    init_states=y_mechanics,
    parameters=p_mechanics,
    num_states=len(y_mechanics),
    fun_monitor=mon_mechanics,
    Ta_index=mechanics_model["monitor_index"]("Ta"),
    missing_variables=missing_mech.values_mechanics,
    num_missing_variables=len(mechanics_ode.missing_variables),
)

problem = mechanicssolver.MechanicsProblem(
    geometry, material, bcs, ode=mechanics_odesolver
)
problem.solve(t0=0.0, dt=0.01)

# activation.vector()[:] = 0.05299967
# problem.solve(t0=0.0, dt=0.01)

# breakpoint()


inds = []
j = 0
for i, ti in enumerate(t):
    print(f"Solving time {ti:.2f} ms")
    ep_solver.step((ti, ti + dt))

    # Update missing values for the mechanics model
    missing_mech.values_ep[0, :] = mv_ep(
        ti, ode._values, ode.parameters, missing_ep.values_ep
    )
    print(f"Missing values (mechanics - ep mesh): {missing_mech.values_ep}")
    missing_mech.ep_values_to_function()
    missing_mech.interpolate_ep_to_mechanics()
    missing_mech.mechanics_function_to_values()
    print(
        f"Missing values (mechanics - mechanics mesh): {missing_mech.values_mechanics}"
    )

    # print(f"Missing values (mechanics): {mechanics_ep_missing_values}")
    # print(f"Missing values (ep): {ep_ep_missing_values}")

    # missing_ep_mech.vector().set_local(mechanics_ep_missing_values[0])
    # missing_ep_mech.vector().set_local(ep_mech_missing_values[0])

    # missing_mechanics_mech.interpolate(missing_ep_mech)
    # mechanics_mech_missing_values[:] = missing_mechanics_mech.vector().get_local()
    # mechanics_ep_missing_values[:] = missing_mechanics_ep.vector().get_local()

    # Compute the change in the missing values
    # change = np.linalg.norm(
    #     mechanics_mech_missing_values - prev_mechanics_missing_values
    # ) / np.linalg.norm(prev_mechanics_missing_values)

    # Check if the change is small enough to continue to the next time step
    # if change < tol:
    #     # Very small change to just continue to next time step
    #     continue
    # breakpoint()
    # Store the index of the time step where we performed a step
    inds.append(i)

    problem.solve(ti, dt)

    if i % 10 == 0:
        U, p = problem.state.split(deepcopy=True)
        with dolfin.XDMFFile("disp.xdmf") as file:
            file.write_checkpoint(U, "disp", j, dolfin.XDMFFile.Encoding.HDF5, True)
        j += 1

    # Forward step for the mechanics model
    # y_mechanics[:] = fgr_mechanics(y_mechanics, ti, dt, p_mechanics, mechanics_missing_values)
    monitor_mechanics = mon_mechanics(
        ti,
        problem.ode._values,
        problem.ode.parameters,
        missing_mech.values_mechanics,
    )

    print(
        f"Ta (mechanics): {monitor_mechanics[mechanics_model['monitor_index']('Ta')]}"
    )
    if i > 100:
        breakpoint()
    # J_TRPN_mechanics[i] = monitor_mechanics[J_TRPN_index_mechanics]

    # Update missing values for the EP model
    missing_ep.values_mechanics[0, :] = mv_mechanics(
        ti,
        problem.ode._values,
        problem.ode.parameters,
        missing_mech.values_mechanics,
    )
    print(f"Missing values (ep - mechanics mesh): {missing_ep.values_mechanics}")
    missing_ep.mechanics_values_to_function()
    missing_ep.interpolate_mechanics_to_ep()
    missing_ep.ep_function_to_values()
    print(f"Missing values (ep - ep mesh): {missing_ep.values_ep}")
    # missing_mechanics_ep.vector().set_local(ep_mech_missing_values[0])
    # missing_ep_ep.interpolate(missing_mechanics_ep)
    # ep_ep_missing_values[:] = missing_ep_ep.vector().get_local()

    # prev_mechanics_missing_values[:] = mechanics_mech_missing_values
    # if i > 500:
    #     break
breakpoint()
# Plot the results
print(f"Solved on {100 * len(inds) / len(t)}% of the time steps")
inds = np.array(inds)

ax[0, 0].plot(t, V_ep, label=f"tol={tol}")
ax[1, 0].plot(t[inds], Ta_mechanics[inds], label=f"tol={tol}")

ax[0, 1].plot(t, Ca_ep, label=f"tol={tol}")

ax[1, 1].plot(t[inds], J_TRPN_mechanics[inds], label=f"tol={tol}")

err_Ta = np.linalg.norm(Ta_full[inds] - Ta_mechanics[inds]) / np.linalg.norm(
    Ta_mechanics
)
err_J_TRPN = np.linalg.norm(
    J_TRPN_full[inds] - J_TRPN_mechanics[inds]
) / np.linalg.norm(J_TRPN_mechanics)
ax[2, 0].plot(
    t[inds],
    Ta_full[inds] - Ta_mechanics[inds],
    label=f"err={err_Ta:.2e}, tol={tol}",
)

ax[2, 1].plot(
    t[inds],
    J_TRPN_full[inds] - J_TRPN_mechanics[inds],
    label=f"err={err_J_TRPN:.2e}, tol={tol}",
)

ax[1, 0].set_xlabel("Time (ms)")
ax[1, 1].set_xlabel("Time (ms)")
ax[0, 0].set_ylabel("V (mV)")
ax[1, 0].set_ylabel("Ta (kPa)")
ax[0, 1].set_ylabel("Ca (mM)")
ax[1, 1].set_ylabel("J TRPN (mM)")

ax[2, 0].set_ylabel("Ta error (kPa)")
ax[2, 1].set_ylabel("J TRPN error (mM)")

for axi in ax.flatten():
    axi.legend()

fig.tight_layout()
fig.savefig("V_and_Ta.png")
