# Heavily ChatGPT-made plots as mentioned in report.
import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ---- Matplotlib defaults ----
# Adjusted for visibility in reports and due to triple figures in report
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 15,
    'figure.figsize': (6, 4),
    'axes.titlesize': 17,
    'axes.labelsize': 20,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.0,
    'legend.fontsize': 17,
    'figure.dpi': 300,
})


os.makedirs("data/plot", exist_ok=True)

# ----- Physical parameters ------
q = 1.0
m = 40.078 # Ca+
z0 = 20.0
d = 500.0
V0 = 25.0e-3 * 9.64852558e7
B = 9.64852558e1
omega_z = np.sqrt(2.0 * q * V0 / (m * d**2))
omega_0 = q * B / m
DEN_TOL = 1e-10  # avoid division by small analytic values in relative errors


def f(t, return_R=False):
    # initial conditions for radial motion
    x0 = 20.0
    v_0y = 25.0

    disc = omega_0**2 - 2.0 * omega_z**2
    sqrt_disc = np.sqrt(disc)
    omega_plus = 0.5 * (omega_0 + sqrt_disc)
    omega_minus = 0.5 * (omega_0 - sqrt_disc)

    print(omega_plus/omega_minus)

    # complex-amplitude solution for f(t) = x + i y
    A_plus  = (v_0y + omega_minus * x0) / (omega_minus - omega_plus)
    A_minus = -(v_0y + omega_plus  * x0) / (omega_minus - omega_plus)
    print(f"A_plus: {A_plus}, A_minus: {A_minus}")
    R_plus = np.abs(np.abs(A_plus) + np.abs(A_minus))
    R_minus = np.abs(np.abs(A_plus) - np.abs(A_minus))
    print(f"R_plus: {R_plus}, R_minus: {R_minus}")
    f_t = A_plus * np.exp(-1j * (omega_plus  * t)) + \
          A_minus* np.exp(-1j * (omega_minus * t))

    x = f_t.real
    y = f_t.imag
    if return_R:
        return x, y, R_plus, R_minus
    else:
        return x, y

# ---- I/O helpers ----
def load_posvel(coulomb_flag: int, particle_idx: int, N: int = 100000): # assumes recommended naming
    path = f"data/pos_vel_{particle_idx}_coulomb={coulomb_flag}_N{N}.txt"
    arr = np.loadtxt(path)
    t, x, y, z, vx, vy, vz = (arr[:,0], arr[:,1], arr[:,2], arr[:,3],
                              arr[:,4], arr[:,5], arr[:,6])
    return dict(t=t, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, path=path)

# Load both particles, with and without Coulomb
p1_on  = load_posvel(1, 0) # coulomb on, particle 0
p2_on  = load_posvel(1, 1)
p1_off = load_posvel(0, 0)
p2_off = load_posvel(0, 1)

print(p1_off["t"])


# ---- Plotting helpers ----
def mark_start_end(x, y, label_prefix=None):
    """Add start (o) and end (x) markers to the current axes."""
    lbl_start = None if label_prefix is None else f"{label_prefix} start"
    lbl_end   = None if label_prefix is None else f"{label_prefix} end"
    plt.plot(x[0],  y[0],  marker="o", markersize=6, linestyle="None", label=lbl_start)
    plt.plot(x[-1], y[-1], marker="x", markersize=6, linestyle="None", label=lbl_end)

def plot_time_series(t1, a1, t2, a2, comp, ylabel, outname, title_suffix="(Coulomb on)"):
    plt.figure()
    plt.plot(t1, a1, alpha=0.9, label=f"P1 {comp}")
    plt.plot(t2, a2, alpha=0.9, label=f"P2 {comp}")
    mark_start_end(t1, a1, "P1")
    mark_start_end(t2, a2, "P2")
    plt.xlabel(r"$t~(\text{µ}\mathrm{s})$")
    plt.ylabel(ylabel)
    # plt.title(f"{comp} vs time {title_suffix}")
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()

def plot_phase_space(x1, v1, x2, v2, xlabel, ylabel, title, outname):
    plt.figure()
    plt.plot(x1, v1, alpha=0.9, label="P1")
    plt.plot(x2, v2, alpha=0.9, label="P2")
    mark_start_end(x1, v1, "P1")
    mark_start_end(x2, v2, "P2")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.title(title)
    plt.axis("equal")
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()

print((p1_off["x"]),(p1_off["vx"]))

def plot_xy(p1, p2, title, outname, legend_on=True):
    """Plot x–y trajectory for both particles."""
    plt.figure()
    plt.plot(p1["x"], p1["y"], alpha=0.9, label="P1")
    plt.plot(p2["x"], p2["y"], alpha=0.9, label="P2")

    # plot circle at R_plus and R_minus
    _, _, R_plus, R_minus = f(0, return_R=True)
    theta = np.linspace(0, 2 * np.pi, 100)
    x_circle_plus = R_plus * np.cos(theta)
    y_circle_plus = R_plus * np.sin(theta)
    plt.plot(x_circle_plus, y_circle_plus, linestyle="--", color="gray", alpha=0.4, label=r"$R_{+}$")
    x_circle_minus = R_minus * np.cos(theta)
    y_circle_minus = R_minus * np.sin(theta)
    plt.plot(x_circle_minus, y_circle_minus, linestyle="--", color="black", alpha=0.4, label=r"$R_{-}$")

    mark_start_end(p1["x"], p1["y"], "P1")
    mark_start_end(p2["x"], p2["y"], "P2")
    plt.xlabel(r"$x~(\text{µ}\mathrm{m})$")
    plt.ylabel(r"$y~(\text{µ}\mathrm{m})$")
    # plt.title(title)
    plt.gca().set_aspect('equal', adjustable='box') # plt.axis("equal") did not work
    if legend_on:
        plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()


# ---- Time-series (Coulomb ON only; names clearer) ----
plot_time_series(p1_on["t"], p1_on["x"], p2_on["t"], p2_on["x"],
                 r"$x$-position", r"$x~(\text{µ}\mathrm{m})$",
                 "data/plot/two_particles_x_vs_time_coulomb_on.pdf")

plot_time_series(p1_on["t"], p1_on["y"], p2_on["t"], p2_on["y"],
                 r"$y$-position", r"$y~(\text{µ}\mathrm{m})$",
                 "data/plot/two_particles_y_vs_time_coulomb_on.pdf")

plot_time_series(p1_on["t"], p1_on["z"], p2_on["t"], p2_on["z"],
                 r"$z$-position", r"$z~(\text{µ}\mathrm{m})$",
                 "data/plot/two_particles_z_vs_time_coulomb_on.pdf")

# ---- Phase space: x–vx and z–vz, with and without Coulomb ----
# Coulomb ON
plot_phase_space(p1_on["x"], p1_on["vx"], p2_on["x"], p2_on["vx"],
                 r"$x~(\text{µ}\mathrm{m})$", r"$v_x~(\text{µ}\mathrm{m}/\text{µ}\mathrm{s})$",
                 r"Phase space for $x$ (Coulomb ON)",
                 "data/plot/phase_x_vx_coulomb_on.pdf")

plot_phase_space(p1_on["z"], p1_on["vz"], p2_on["z"], p2_on["vz"],
                 r"$z~(\text{µ}\mathrm{m})$", r"$v_z~(\text{µ}\mathrm{m}/\text{µ}\mathrm{s})$",
                 r"Phase space for $z$ (Coulomb ON)",
                 "data/plot/phase_z_vz_coulomb_on.pdf")

# Coulomb OFF
plot_phase_space(p1_off["x"], p1_off["vx"], p2_off["x"], p2_off["vx"],
                 r"$x~(\text{µ}\mathrm{m})$", r"$v_x~(\text{µ}\mathrm{m}/\text{µ}\mathrm{s})$",
                 r"Phase space for $x$ (Coulomb OFF)",
                 "data/plot/phase_x_vx_coulomb_off.pdf")

plot_phase_space(p1_off["z"], p1_off["vz"], p2_off["z"], p2_off["vz"],
                 r"$z~(\text{µ}\mathrm{m})$", r"$v_z~(\text{µ}\mathrm{m}/\text{µ}\mathrm{s})$",
                 r"Phase space for $z$ (Coulomb OFF)",
                 "data/plot/phase_z_vz_coulomb_off.pdf")



# ---- XY Trajectories ----
plot_xy(p1_on, p2_on, r"XY trajectory (Coulomb ON)",  "data/plot/xy_traj_coulomb_on.pdf")
plot_xy(p1_off, p2_off, r"XY trajectory (Coulomb OFF)", "data/plot/xy_traj_coulomb_off.pdf")

# plot cover page
p1_on_cover = load_posvel(0, 0, N=1000000) # coulomb on, particle 0, longer run
p2_on_cover = load_posvel(0, 1, N=1000000) # coulomb on, particle 1, longer run

plt.figure()
plt.plot(p1_on_cover["x"], p1_on_cover["y"], alpha=0.9, label="P1")
plt.plot(p2_on_cover["x"], p2_on_cover["y"], alpha=0.9, label="P2")

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

plt.tight_layout()
plt.savefig("data/plot/xy_traj_coulomb_off_cover.pdf", dpi=300, bbox_inches="tight")
plt.close()

# Plot energy
energy_data0 = np.loadtxt("data/total_energy_coulomb=0_N40000.txt")
energy_data1 = np.loadtxt("data/total_energy_coulomb=1_N40000.txt")
t  = energy_data0[:, 0]
e0 = energy_data0[:, 1]
e1 = energy_data1[:, 1]

E0_0 = e0[0]
E0_1 = e1[0]
fig, ax = plt.subplots()
ax.plot(t, e1 - E0_1, label="Coulomb ON")
ax.plot(t, e0 - E0_0, label="Coulomb OFF")

ax.set_xlabel(r"$t\;(\mu\mathrm{s})$")
ax.set_ylabel(r"$\Delta E\;[u\,(m/s)^2]$")
ax.axhline(0, lw=1, alpha=0.6, color="gray", linestyle="--", label="Initial Energy")

ax.ticklabel_format(axis='y', style='plain', useOffset=False)
ax.yaxis.get_offset_text().set_visible(False)
ax.legend(frameon=True)
ax.margins(y=0.05)
plt.grid(True, which="both", alpha=0.3, linewidth=0.6, linestyle="--")
plt.tight_layout()
plt.savefig("data/plot/total_energy_delta_few_particles.pdf")
plt.close()

# Plot relative error in energy
rel_err0 = np.abs(e0 - e0[0]) / (np.abs(e0[0]) + DEN_TOL)
rel_err1 = np.abs(e1 - e1[0]) / (np.abs(e1[0]) + DEN_TOL)
plt.figure()
plt.plot(t, rel_err1, label="Coulomb ON")
plt.plot(t, rel_err0, label="Coulomb OFF")
plt.xlabel(r"$t~(\text{µ}\mathrm{s})$")
plt.ylabel(r"Relative error in Total Energy")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig("data/plot/total_energy_relerr_few_particles.pdf")