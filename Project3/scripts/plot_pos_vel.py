# Heavily ChatGPT-made plots as mentioned in report.
import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

# ---- Matplotlib defaults ----
# Adjusted for visibility in reports and due to triple figures in report
plt.rcParams.update({
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

def plot_xy(p1, p2, title, outname):
    """Plot x–y trajectory for both particles."""
    plt.figure()
    plt.plot(p1["x"], p1["y"], alpha=0.9, label="P1")
    plt.plot(p2["x"], p2["y"], alpha=0.9, label="P2")
    mark_start_end(p1["x"], p1["y"], "P1")
    mark_start_end(p2["x"], p2["y"], "P2")
    plt.xlabel(r"$x~(\text{µ}\mathrm{m})$")
    plt.ylabel(r"$y~(\text{µ}\mathrm{m})$")
    # plt.title(title)
    plt.gca().set_aspect('equal', adjustable='box') # plt.axis("equal") did not work
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