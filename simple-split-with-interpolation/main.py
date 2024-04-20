"""Coupled model without interpolating the variables between
the mechanics and ep mesh. Instead we just use the means value. 
This mean that we cannot really model any spatial variation.
"""

import gotranx
from pathlib import Path
from typing import Any
import numpy as np
import dolfin
import pulse
import beat
import ufl_legacy as ufl
import matplotlib.pyplot as plt
import mechanicssolver
import interpolation

try:
    raise ImportError
    from numba import jit
except ImportError:

    def jit(*args, **kwargs):
        def wrapper(func):
            return func

        return wrapper


dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["quadrature_degree"] = 4


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
if not Path("ep_model.py").exists():
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

    Path("ep_model.py").write_text(code_ep)
    Path("mechanics_model.py").write_text(code_mechanics)

import ep_model as _ep_model
import mechanics_model as _mechanics_model

ep_model = _ep_model.__dict__
mechanics_model = _mechanics_model.__dict__

# Set time step to 0.1 ms
dt = 0.05
# Simulate model for 1000 ms
t = np.arange(0, 200, dt)

# Get the index of the membrane potential
V_index_ep = ep_model["state_index"]("v")
# Forwared generalized rush larsen scheme for the electrophysiology model
fgr_ep = jit(nopython=True)(ep_model["forward_generalized_rush_larsen"])
# Monitor function for the electrophysiology model
mon_ep = ep_model["monitor_values"]
# Missing values function for the electrophysiology model
mv_ep = ep_model["missing_values"]
# Index of the calcium concentration
Ca_index_ep = ep_model["state_index"]("cai")

# Forwared generalized rush larsen scheme for the mechanics model
fgr_mechanics = jit(nopython=True)(mechanics_model["forward_generalized_rush_larsen"])
# Monitor function for the mechanics model
mon_mechanics = mechanics_model["monitor_values"]
# Missing values function for the mechanics model
mv_mechanics = mechanics_model["missing_values"]
# Index of the active tension
Ta_index_mechanics = mechanics_model["monitor_index"]("Ta")

# lmbda_index = mechanics_model["parameter_index"]("lambda")
# dLmbda_idnex = mechanics_model["parameter_index"]("lambda_rate")

lmbda_index = mechanics_model["parameter_index"]("lmbda")
dLmbda_idnex = mechanics_model["parameter_index"]("dLambda")
lmbda_index_ep = ep_model["parameter_index"]("lmbda")
# Index of the J_TRPN
# J_TRPN_index_mechanics = mechanics_model["monitor_index"]("J_TRPN")


tol = 5e-4
# Create arrays to store the results
V_ep = np.zeros(len(t))
Ca_ep = np.zeros(len(t))

Ta_mechanics = np.zeros(len(t))
J_TRPN_mechanics = np.zeros(len(t))


# Get initial values from the EP model
y_ep_ = ep_model["init_state_values"]()
p_ep_ = ep_model["init_parameter_values"](amp=0.0)

ep_missing_values_ = np.zeros(len(ep_model["missing"]))

# Get initial values from the mechanics model
y_mechanics_ = mechanics_model["init_state_values"]()
p_mechanics_ = mechanics_model["init_parameter_values"]()
mechanics_missing_values_ = np.zeros(len(mechanics_model["missing"]))

mesh = setup_geometry(dx=3.0)
ep_mesh = dolfin.adapt(dolfin.adapt(dolfin.adapt(mesh)))


# Surface to volume ratio
chi = 140.0  # mm^{-1}
# Membrane capacitance
C_m = 0.01  # mu F / mm^2

time = dolfin.Constant(0.0)
I_s = define_stimulus(mesh=ep_mesh, chi=chi, C_m=C_m, time=time)

M = define_conductivity_tensor(chi, C_m)
# M = ufl.zero((3, 3))

# params = {"linear_solver_type": "direct"}
params = {"preconditioner": "sor", "use_custom_preconditioner": False}
ep_ode_space = dolfin.FunctionSpace(ep_mesh, "CG", 1)
v_ode = dolfin.Function(ep_ode_space)
num_points_ep = v_ode.vector().local_size()
lmbda = dolfin.Function(ep_ode_space)

y_ep = np.zeros((len(y_ep_), num_points_ep))
y_ep.T[:] = y_ep_


# Get the default values of the missing values
# A little bit chicken and egg problem here, but in this specific case we know that
# the mechanics_missing_values is only the calcium concentration, which is a state variable
# and this doesn't require any additional information to be calculated.
mechanics_missing_values_[:] = mv_ep(0, y_ep_, p_ep_, ep_missing_values_)
ep_missing_values_[:] = mv_mechanics(
    0, y_mechanics_, p_mechanics_, mechanics_missing_values_
)


# Set the activation
# activation_space = dolfin.FunctionSpace(mesh, "CG", 1)
activation_space = dolfin.FunctionSpace(mesh, "CG", 1)
activation = dolfin.Function(activation_space)
num_points_mech = activation.vector().local_size()


missing_mech = interpolation.MissingValue(
    element=activation.ufl_element(),
    interpolation_element=ep_ode_space.ufl_element(),
    mechanics_mesh=mesh,
    ep_mesh=ep_mesh,
    num_values=len(mechanics_missing_values_),
)

missing_ep = interpolation.MissingValue(
    element=ep_ode_space.ufl_element(),
    interpolation_element=activation.ufl_element(),
    mechanics_mesh=mesh,
    ep_mesh=ep_mesh,
    num_values=len(ep_missing_values_),
)


missing_mech.values_ep.T[:] = mechanics_missing_values_
missing_ep.values_mechanics.T[:] = ep_missing_values_
missing_ep.values_ep.T[:] = ep_missing_values_
missing_mech.values_mechanics.T[:] = mechanics_missing_values_

p_ep = np.zeros((len(p_ep_), num_points_ep))
p_ep.T[:] = p_ep_


pde = beat.MonodomainModel(time=time, mesh=ep_mesh, M=M, I_s=I_s, params=params)
ode = beat.odesolver.DolfinODESolver(
    v_ode=dolfin.Function(ep_ode_space),
    v_pde=pde.state,
    fun=fgr_ep,
    parameters=p_ep,
    init_states=y_ep,
    num_states=len(y_ep),
    v_index=ep_model["state_index"]("v"),
    missing_variables=missing_ep.values_ep,
    num_missing_variables=len(ep_model["missing"]),
)

ep_solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode, theta=0.5)

V_ep = np.zeros(len(t))
Ca_ep = np.zeros(len(t))
Ta_mechanics = np.zeros(len(t))
J_TRPN_mechanics = np.zeros(len(t))


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
# material_parameters = pulse.HolzapfelOgden.default_parameters()

# Select model for active contraction
# active_model = "active_stress"


y_mechanics = np.zeros((len(y_mechanics_), num_points_mech))
y_mechanics.T[:] = y_mechanics_
p_mechanics = np.zeros((len(p_mechanics_), num_points_mech))
p_mechanics.T[:] = p_mechanics_


mechanics_missing_values = np.zeros(
    (len(mechanics_model["missing"]), activation.vector().local_size())
)
mechanics_missing_values.T[:] = mechanics_missing_values_


material_parameters = dict(
    a=2.28,
    a_f=1.686,
    b=9.726,
    b_f=15.779,
    a_s=0.0,
    b_s=0.0,
    a_fs=0.0,
    b_fs=0.0,
)
material = pulse.HolzapfelOgden(
    active_model="active_stress",
    parameters=material_parameters,
    activation=activation,
    f0=f0,
    s0=s0,
    n0=n0,
)


# Make Dirichlet boundary conditions
def dirichlet_bc(W):
    V = W if W.sub(0).num_sub_spaces() == 0 else W.sub(0)
    return dolfin.DirichletBC(V, dolfin.Constant((0.0, 0.0, 0.0)), fixed)


# Collect Boundary Conditions
bcs = pulse.BoundaryConditions(dirichlet=(dirichlet_bc,))

mech_ode = mechanicssolver.ODESolver(
    Ta=activation,
    fun=fgr_mechanics,
    states=y_mechanics,
    parameters=p_mechanics,
    fun_monitor=mon_mechanics,
    Ta_index=Ta_index_mechanics,
    missing_variables=missing_mech.values_mechanics,
    num_missing_variables=len(mechanics_missing_values),
    lmbda_index=lmbda_index,
    dLambda_index=dLmbda_idnex,
    f0=f0,
)

problem = mechanicssolver.MechanicsProblem(geometry, material, bcs, ode=mech_ode)
problem.solve(0.0, 0.0)


disp_file = Path("disp.xdmf")
disp_file.unlink(missing_ok=True)
disp_file.with_suffix(".h5").unlink(missing_ok=True)

V_file = Path("V.xdmf")
V_file.unlink(missing_ok=True)
V_file.with_suffix(".h5").unlink(missing_ok=True)


Ta_file = Path("Ta.xdmf")
Ta_file.unlink(missing_ok=True)
Ta_file.with_suffix(".h5").unlink(missing_ok=True)

lmbda_file = Path("lmbda.xdmf")
lmbda_file.unlink(missing_ok=True)
lmbda_file.with_suffix(".h5").unlink(missing_ok=True)
import logging


logging.getLogger("beat").setLevel(logging.ERROR)

inds = []
j = 0
N = 10
theta = 0.5
for i, ti in enumerate(t):
    print(f"Solving time {ti:.2f} ms")
    ep_solver.step((ti, ti + dt))

    missing_ep_values = mv_ep(
        ti + dt, ode._values, ode.parameters, missing_ep.values_ep
    )

    for k in range(missing_mech.num_values):
        missing_mech.u_ep_int[k].vector()[:] = missing_ep_values[k, :]
    # missing_mech.interpolate_ep()

    # print(f"Missing values (mechanics - ep mesh): {missing_mech.values_ep}")
    # missing_mech.ep_values_to_function()
    missing_mech.interpolate_ep_to_mechanics()
    missing_mech.mechanics_function_to_values()
    print(
        f"Missing values (mechanics - mechanics mesh): min: {missing_mech.values_mechanics.min()}, max: {missing_mech.values_mechanics.max()}"
    )

    # breakpoint()

    vi = ode._values[V_index_ep]
    # print(f"V: {vi[0]}, min(V): {vi.min()}, max(V): {vi.max()}")
    V_ep[i] = ode._values[V_index_ep][0]
    Ca_ep[i] = ode._values[Ca_index_ep][0]

    if i % N != 0:
        continue
    inds.append(i)

    # mech_ode.step(ti, dt)
    problem.solve(ti, dt * N)

    missing_mech_values = mv_mechanics(
        ti + dt * N,
        mech_ode.states,
        mech_ode.parameters,
        missing_mech.values_mechanics,
    )

    for k in range(missing_ep.num_values):
        missing_ep.u_mechanics_int[k].vector()[:] = missing_mech_values[k, :]

    # missing_ep.interpolate_mechanics()
    # print(f"Missing values (ep - mechanics mesh): {missing_ep.values_mechanics}")
    # missing_ep.mechanics_values_to_function()
    missing_ep.interpolate_mechanics_to_ep()
    missing_ep.ep_function_to_values()
    lmbda.interpolate(problem.solver.lmbda)
    p_ep[lmbda_index_ep, :] = lmbda.vector().get_local()

    # print(f"Missing values (ep - ep mesh): {missing_ep.values_ep}")

    # if i % 10 == 0:
    U, p = problem.state.split(deepcopy=True)
    with dolfin.XDMFFile(disp_file.as_posix()) as file:
        file.write_checkpoint(U, "disp", j, dolfin.XDMFFile.Encoding.HDF5, True)
    with dolfin.XDMFFile(V_file.as_posix()) as file:
        file.write_checkpoint(pde.state, "V", j, dolfin.XDMFFile.Encoding.HDF5, True)
    with dolfin.XDMFFile(Ta_file.as_posix()) as file:
        file.write_checkpoint(activation, "Ta", j, dolfin.XDMFFile.Encoding.HDF5, True)
    with dolfin.XDMFFile(lmbda_file.as_posix()) as file:
        file.write_checkpoint(
            problem.solver.lmbda, "lambda", j, dolfin.XDMFFile.Encoding.HDF5, True
        )
    j += 1


# Plot the results
print(f"Solved on {100 * len(inds) / len(t)}% of the time steps")
inds = np.array(inds)
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(10, 10))
ax[0, 0].plot(t, V_ep, label=f"tol={tol}")
ax[1, 0].plot(t[inds], Ta_mechanics[inds], label=f"tol={tol}")

ax[0, 1].plot(t, Ca_ep, label=f"tol={tol}")

ax[1, 1].plot(t[inds], J_TRPN_mechanics[inds], label=f"tol={tol}")


ax[1, 0].set_xlabel("Time (ms)")
ax[1, 1].set_xlabel("Time (ms)")
ax[0, 0].set_ylabel("V (mV)")
ax[1, 0].set_ylabel("Ta (kPa)")
ax[0, 1].set_ylabel("Ca (mM)")
ax[1, 1].set_ylabel("J TRPN (mM)")


for axi in ax.flatten():
    axi.legend()

fig.tight_layout()
fig.savefig("V_and_Ta.png")
