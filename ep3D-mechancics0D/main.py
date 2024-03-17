"""Solve 3D EP model and 0D mechanics model"""

import gotranx
from typing import Any
import numpy as np
import dolfin
import beat
import ufl_legacy as ufl
import matplotlib.pyplot as plt


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
t = np.arange(0, 50, dt)

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


# Get initial values from the EP model
y_ep_ = ep_model["init_state_values"]()
p_ep = ep_model["init_parameter_values"]()  # amp=0.0)
ep_missing_values_ = np.zeros(len(ep_ode.missing_variables))

# Get initial values from the mechanics model
y_mechanics = mechanics_model["init_state_values"]()
p_mechanics = mechanics_model["init_parameter_values"]()
mechanics_missing_values = np.zeros(len(mechanics_ode.missing_variables))


mesh = setup_geometry(dx=0.5)

# Surface to volume ratio
chi = 140.0  # mm^{-1}
# Membrane capacitance
C_m = 0.01  # mu F / mm^2

time = dolfin.Constant(0.0)
I_s = define_stimulus(mesh=mesh, chi=chi, C_m=C_m, time=time, A=0)

# M = define_conductivity_tensor(chi, C_m)
M = ufl.zero((3, 3))

# params = {"linear_solver_type": "direct"}
params = {"preconditioner": "sor", "use_custom_preconditioner": False}
ep_ode_space = dolfin.FunctionSpace(mesh, "Lagrange", 1)
v_ode = dolfin.Function(ep_ode_space)
num_points_ep = v_ode.vector().local_size()

y_ep = np.zeros((len(y_ep_), num_points_ep))
y_ep.T[:] = y_ep_


# Get the default values of the missing values
# A little bit chicken and egg problem here, but in this specific case we know that
# the mechanics_missing_values is only the calcium concentration, which is a state variable
# and this doesn't require any additional information to be calculated.
mechanics_missing_values[:] = mv_ep(0, y_ep_, p_ep, ep_missing_values_)
ep_missing_values_[:] = mv_mechanics(
    0, y_mechanics, p_mechanics, mechanics_missing_values
)

ep_missing_values = np.zeros((len(ep_ode.missing_variables), num_points_ep))
ep_missing_values.T[:] = ep_missing_values_


pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
ode = beat.odesolver.DolfinODESolver(
    v_ode=dolfin.Function(ep_ode_space),
    v_pde=pde.state,
    fun=fgr_ep,
    init_states=y_ep,
    parameters=p_ep,
    num_states=len(y_ep),
    v_index=ep_model["state_index"]("v"),
    missing_variables=ep_missing_values,
    num_missing_variables=len(ep_ode.missing_variables),
)

ep_solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode, theta=0.5)

V_ep = np.zeros(len(t))
Ca_ep = np.zeros(len(t))
Ta_mechanics = np.zeros(len(t))
J_TRPN_mechanics = np.zeros(len(t))


inds = []
j = 0
theta = 0.5
for i, ti in enumerate(t):
    print(f"Solving time {ti:.2f} ms")
    ep_solver.step((ti, ti + dt))

    # Update missing values for the mechanics model
    mechanics_missing_values_ep = mv_ep(
        ti, ode._values, ode.parameters, ep_missing_values
    )
    mechanics_missing_values[:] = mechanics_missing_values_ep.mean()

    # Just pick the first value
    vi = ode._values[V_index_ep]
    print(f"V: {vi[0]}, min(V): {vi.min()}, max(V): {vi.max()}")
    V_ep[i] = ode._values[V_index_ep][0]
    Ca_ep[i] = ode._values[Ca_index_ep][0]

    inds.append(i)

    y_mechanics[:] = fgr_mechanics(
        y_mechanics, ti, dt, p_mechanics, mechanics_missing_values
    )

    monitor_mechanics = mon_mechanics(
        ti,
        y_mechanics,
        p_mechanics,
        mechanics_missing_values,
    )
    Ta_mechanics[i] = monitor_mechanics[Ta_index_mechanics]
    J_TRPN_mechanics[i] = monitor_mechanics[J_TRPN_index_mechanics]

    ep_missing_values.T[:] = mv_mechanics(
        ti, y_mechanics, p_mechanics, mechanics_missing_values
    )

    print(f"Ta (mechanics): {Ta_mechanics[i]}")


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
