"""Solve all ODEs on the EP side and simply transfer the activation to the mechanics model"""

import gotranx
from typing import Any
import numpy as np
import dolfin
import beat
import pulse
from pathlib import Path


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


def define_stimulus(mesh, chi, C_m, time):
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
    duration = 2.0  # ms
    A = 50000.0  # mu A/cm^3
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

# Generate code for the electrophysiology model
code = gotranx.cli.gotran2py.get_code(
    ode,
    scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen],
)

# # Execute code and get the models
model: dict[str, Any] = {}
exec(code, model)


# Set time step to 0.1 ms
dt = 0.1
# Simulate model for 1000 ms
t = np.arange(0, 500, dt)

# Get the index of the membrane potential
V_index_ep = model["state_index"]("v")
# Forwared generalized rush larsen scheme for the electrophysiology model
fgr = model["forward_generalized_rush_larsen"]
# Monitor function for the electrophysiology model
mon = model["monitor"]
# Index of the calcium concentration
Ca_index = model["state_index"]("cai")
# Index of the active tension
Ta_index = model["monitor_index"]("Ta")

tol = 2e-5
# Create arrays to store the results
V_ep = np.zeros(len(t))
Ca_ep = np.zeros(len(t))

Ta_mechanics = np.zeros(len(t))

# Get initial values from the EP model
init_states = model["init_state_values"]()
parameters = model["init_parameter_values"](amp=0.0)

mesh = setup_geometry(dx=3.0)
ep_mesh = dolfin.adapt(dolfin.adapt(dolfin.adapt(mesh)))

# Surface to volume ratio
chi = 140.0  # mm^{-1}
# Membrane capacitance
C_m = 0.01  # mu F / mm^2

time = dolfin.Constant(0.0)
I_s = define_stimulus(mesh=ep_mesh, chi=chi, C_m=C_m, time=time)

M = define_conductivity_tensor(chi, C_m)

# params = {"linear_solver_type": "direct"}
params = {"preconditioner": "sor", "use_custom_preconditioner": False}
ode_space = dolfin.FunctionSpace(ep_mesh, "Lagrange", 1)
v_ode = dolfin.Function(ode_space)
num_points = v_ode.vector().local_size()

y = np.zeros((len(init_states), num_points))
y.T[:] = init_states


pde = beat.MonodomainModel(time=time, mesh=ep_mesh, M=M, I_s=I_s, params=params)
ode = beat.odesolver.DolfinODESolver(
    v_ode=dolfin.Function(ode_space),
    v_pde=pde.state,
    fun=fgr,
    init_states=y,
    parameters=parameters,
    num_states=len(y),
    v_index=model["state_index"]("v"),
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
# active_model = pulse.ActiveModels.active_stress
active_model = "active_stress"

# Set the activation
V = dolfin.FunctionSpace(geometry.mesh, "Lagrange", 1)
activation = dolfin.Function(V)
V_ep = dolfin.FunctionSpace(ep_mesh, "Lagrange", 1)
activation_ep = dolfin.Function(V_ep)

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


problem = pulse.MechanicsProblem(geometry, material, bcs)
problem.solve()

result_file = Path("results.xdmf")
result_file.unlink(missing_ok=True)
result_file.with_suffix(".h5").unlink(missing_ok=True)
inds = []
j = 0
for i, ti in enumerate(t):
    print(f"Solving time {ti:.2f} ms")
    ep_solver.step((ti, ti + dt))
    monitor = mon(ti, ode._values, ode.parameters)

    activation_ep.vector().set_local(monitor[Ta_index])
    activation.interpolate(activation_ep)
    print(
        "Ta: ",
        activation.vector().get_local().max(),
        activation.vector().get_local().min(),
    )
    print(
        "V: ",
        ep_solver.pde.state.vector().get_local().max(),
        ep_solver.pde.state.vector().get_local().min(),
    )

    if i % 10 == 0:
        problem.solve()
        U, p = problem.state.split(deepcopy=True)
        with dolfin.XDMFFile(str(result_file)) as file:
            file.write_checkpoint(U, "disp", j, dolfin.XDMFFile.Encoding.HDF5, True)
            file.write_checkpoint(
                problem.material.activation,
                "Ta",
                j,
                dolfin.XDMFFile.Encoding.HDF5,
                True,
            )
            file.write_checkpoint(
                pde.state, "V", j, dolfin.XDMFFile.Encoding.HDF5, True
            )
        j += 1
