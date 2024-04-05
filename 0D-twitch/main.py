"""Same as 0D but with varying lambda
"""

import gotranx
from typing import Any
import numpy as np
import matplotlib.pyplot as plt


def twitch(t, tstart=0.05, ca_ampl=-0.2):
    tau1 = 0.05 * 1000
    tau2 = 0.110 * 1000

    ca_diast = 0.0

    beta = (tau1 / tau2) ** (-1 / (tau1 / tau2 - 1)) - (tau1 / tau2) ** (
        -1 / (1 - tau2 / tau1)
    )
    ca = np.zeros_like(t)

    ca[t <= tstart] = ca_diast

    ca[t > tstart] = (ca_ampl - ca_diast) / beta * (
        np.exp(-(t[t > tstart] - tstart) / tau1)
        - np.exp(-(t[t > tstart] - tstart) / tau2)
    ) + ca_diast
    return ca + 1.0


# Load the model
ode = gotranx.load_ode("ORdmm_Land.ode")

mechanics_comp = ode.get_component("mechanics")
mechanics_ode = mechanics_comp.to_ode()

ep_ode = ode - mechanics_comp

# Generate code for full model
code = gotranx.cli.gotran2py.get_code(
    ode,
    scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen],
)

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
model: dict[str, Any] = {}
exec(code, model)
ep_model: dict[str, Any] = {}
exec(code_ep, ep_model)
mechanics_model: dict[str, Any] = {}
exec(code_mechanics, mechanics_model)


# Set time step to 0.1 ms
dt = 0.1
# Simulate model for 1000 ms
t = np.arange(0, 1000, dt)
fig, ax = plt.subplots()
ax.plot(t, twitch(t))
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Lambda")
fig.savefig("twitch.png")

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
lmbda_index_mechanics = mechanics_model["parameter_index"]("lmbda")

# Forwared generalized rush larsen scheme for the full model
fgr = model["forward_generalized_rush_larsen"]
# Monitor function for the full model
mon = model["monitor"]
# Index of the active tension for the full model
Ta_index = model["monitor_index"]("Ta")
# Index of the J_TRPN for the full model
J_TRPN_index = model["monitor_index"]("J_TRPN")
lmbda_index = model["parameter_index"]("lmbda")


# Tolerances to test for when to perform steps in the mechanics model
tols = [
    2e-5,
    4e-5,
    6e-5,
    8e-5,
    2e-4,
    4e-4,
    6e-4,
    8e-4,
    2e-3,
    4e-3,
    6e-3,
    8e-3,
    1e-2,
]
# Colors for the plots
from itertools import cycle

colors = cycle(["r", "g", "b", "c", "m"])
linestyles = cycle(["-", "--", "-.", ":"])

# Create arrays to store the results
V_ep = np.zeros(len(t))
Ca_ep = np.zeros(len(t))

J_TRPN_full = np.zeros(len(t))
Ta_full = np.zeros(len(t))

Ta_mechanics = np.zeros(len(t))
J_TRPN_mechanics = np.zeros(len(t))

fig, ax = plt.subplots(3, 2, sharex=True, figsize=(14, 10))
lines = []
labels = []

for j, (col, ls, tol) in enumerate(zip(colors, linestyles, tols)):
    # Get initial values from the EP model
    y_ep = ep_model["init_state_values"]()
    p_ep = ep_model["init_parameter_values"]()
    ep_missing_values = np.zeros(len(ep_ode.missing_variables))

    # Get initial values from the mechanics model
    y_mechanics = mechanics_model["init_state_values"]()
    p_mechanics = mechanics_model["init_parameter_values"]()
    mechanics_missing_values = np.zeros(len(mechanics_ode.missing_variables))

    # Get the initial values from the full model
    y = model["init_state_values"]()
    p = model["init_parameter_values"]()

    # Get the default values of the missing values
    # A little bit chicken and egg problem here, but in this specific case we know that
    # the mechanics_missing_values is only the calcium concentration, which is a state variable
    # and this doesn't require any additional information to be calculated.
    mechanics_missing_values[:] = mv_ep(0, y_ep, p_ep, ep_missing_values)
    ep_missing_values[:] = mv_mechanics(
        0, y_mechanics, p_mechanics, mechanics_missing_values
    )

    # We will store the previous missing values to check for convergence
    prev_mechanics_missing_values = np.zeros_like(mechanics_missing_values)
    prev_mechanics_missing_values[:] = mechanics_missing_values

    inds = []
    count = 1
    max_count = 10
    for i, ti in enumerate(t):

        # Update lambda
        p[lmbda_index] = twitch(ti)
        p_mechanics[lmbda_index_mechanics] = twitch(ti)

        # Forward step for the full model
        y[:] = fgr(y, ti, dt, p)
        monitor = mon(ti, y, p)
        J_TRPN_full[i] = monitor[J_TRPN_index]
        Ta_full[i] = monitor[Ta_index]

        # Forward step for the EP model
        y_ep[:] = fgr_ep(y_ep, ti, dt, p_ep, ep_missing_values)
        V_ep[i] = y_ep[V_index_ep]
        Ca_ep[i] = y_ep[Ca_index_ep]
        monitor_ep = mon_ep(ti, y_ep, p_ep, ep_missing_values)

        # Update missing values for the mechanics model
        mechanics_missing_values[:] = mv_ep(t, y_ep, p_ep, ep_missing_values)

        # Compute the change in the missing values
        change = np.linalg.norm(
            mechanics_missing_values - prev_mechanics_missing_values
        ) / np.linalg.norm(prev_mechanics_missing_values)

        # Check if the change is small enough to continue to the next time step
        if change < tol:
            count += 1
            # Very small change to just continue to next time step
            if count < max_count:
                continue

        # Store the index of the time step where we performed a step
        inds.append(i)

        # Forward step for the mechanics model
        y_mechanics[:] = fgr_mechanics(
            y_mechanics, ti, count * dt, p_mechanics, mechanics_missing_values
        )
        count = 1
        monitor_mechanics = mon_mechanics(
            ti,
            y_mechanics,
            p_mechanics,
            mechanics_missing_values,
        )
        Ta_mechanics[i] = monitor_mechanics[Ta_index_mechanics]
        J_TRPN_mechanics[i] = monitor_mechanics[J_TRPN_index_mechanics]

        # Update missing values for the EP model
        ep_missing_values[:] = mv_mechanics(
            t, y_mechanics, p_mechanics, mechanics_missing_values
        )

        prev_mechanics_missing_values[:] = mechanics_missing_values

    # Plot the results
    perc = 100 * len(inds) / len(t)
    print(f"Solved on {perc}% of the time steps")
    inds = np.array(inds)

    if j == 0:
        # Plot the full model with a dashed line only for the first run
        (l,) = ax[1, 0].plot(t, Ta_full, color="k", linestyle="--", label="Full")
        ax[1, 1].plot(t, J_TRPN_full, color="k", linestyle="--", label="Full")
        lines.append(l)
        labels.append("Full")

    (l,) = ax[0, 0].plot(t, V_ep, color=col, linestyle=ls, label=f"tol={tol}")
    lines.append(l)
    labels.append(f"tol={tol}, perc={perc}%")
    ax[1, 0].plot(
        t[inds],
        Ta_mechanics[inds],
        color=col,
        linestyle=ls,  # label=f"tol={tol}"
    )

    ax[0, 1].plot(
        t,
        Ca_ep,
        color=col,
        linestyle=ls,
        # label=f"tol={tol}",
    )

    ax[1, 1].plot(
        t[inds],
        J_TRPN_mechanics[inds],
        color=col,
        linestyle=ls,  # label=f"tol={tol}"
    )

    err_Ta = np.linalg.norm(Ta_full[inds] - Ta_mechanics[inds]) / np.linalg.norm(
        Ta_mechanics
    )
    err_J_TRPN = np.linalg.norm(
        J_TRPN_full[inds] - J_TRPN_mechanics[inds]
    ) / np.linalg.norm(J_TRPN_mechanics)

    ax[2, 0].plot(
        t[inds],
        Ta_full[inds] - Ta_mechanics[inds],
        # label=f"tol={tol}, perc={perc}%",
        color=col,
        linestyle=ls,
    )

    ax[2, 1].plot(
        t[inds],
        J_TRPN_full[inds] - J_TRPN_mechanics[inds],
        # label=f"err={err_J_TRPN:.2e}, tol={tol}",
        color=col,
        linestyle=ls,
    )

    ax[1, 0].set_xlabel("Time (ms)")
    ax[1, 1].set_xlabel("Time (ms)")
    ax[0, 0].set_ylabel("V (mV)")
    ax[1, 0].set_ylabel("Ta (kPa)")
    ax[0, 1].set_ylabel("Ca (mM)")
    ax[1, 1].set_ylabel("J TRPN (mM)")

    ax[2, 0].set_ylabel("Ta error (kPa)")
    ax[2, 1].set_ylabel("J TRPN error (mM)")

    # for axi in ax.flatten():
    # axi.legend()

    # fig.tight_layout()
    if j == len(tols) - 1:
        fig.subplots_adjust(right=0.95)
        lgd = fig.legend(lines, labels, loc="center left", bbox_to_anchor=(1.0, 0.5))
        fig.savefig("V_and_Ta.png", bbox_extra_artists=(lgd,), bbox_inches="tight")
    else:
        fig.savefig("V_and_Ta.png")
