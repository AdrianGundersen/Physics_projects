# Heavily chatGPT made plots as mentioned in report.
import os
import numpy as np
import matplotlib.pyplot as plt

# ---- Matplotlib defaults (larger, readable in PDFs) ----
plt.rcParams.update({
    'font.size': 16,
    'figure.figsize': (6, 8),   # default; each figure will be 2 rows stacked
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'lines.linewidth': 2.2,
    'legend.fontsize': 12,
    'figure.dpi': 300,
})

# ---- I/O ----
os.makedirs("data/plot", exist_ok=True)
data = np.loadtxt("data/single_particle.txt")

# Columns: [t  x_rk  y_rk  z_rk  vx_rk  vy_rk  vz_rk  t_eu  x_eu  y_eu  z_eu  vx_eu  vy_eu  vz_eu]
t,  x_rk,  y_rk,  z_rk,  vx_rk,  vy_rk,  vz_rk  = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
t_eu, x_eu, y_eu, z_eu, vx_eu, vy_eu, vz_eu = data[:,7], data[:,8], data[:,9], data[:,10], data[:,11], data[:,12], data[:,13]

# ---- Analytical z(t) ----
q = 1.0
m = 1.0
z0 = 20.0     # µm
d = 500.0     # µm
V0 = 2.41e6   # V

omega_z = np.sqrt(2.0 * q * V0 / (m * d**2))  # rad/µs (given your units)
z_analytic_rk = z0 * np.cos(omega_z * t)
z_analytic_eu = z0 * np.cos(omega_z * t_eu)

def stacked_panels(name: str, t_num, z_num, z_an, outfile: str):
    """
    Top: numerical vs analytical z(t).
    Bottom: residuals (numerical - analytical).
    Stacked vertically for readability in PDFs.
    """
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
    (ax_top, ax_bot) = axes

    # Top: trajectories
    ax_top.plot(t_num, z_num, label=f"{name} z(t)")
    ax_top.plot(t_num, z_an, "--", label="Analytical z(t)")
    ax_top.set_ylabel(r"$z~(\mu \mathrm{m})$")
    ax_top.legend(loc="best")
    ax_top.set_title(f"{name}: z(t) vs analytical")
    ax_top.grid(True, linestyle=":", linewidth=0.8)

    # Bottom: residuals
    res = z_num - z_an
    ax_bot.plot(t_num, res, label=f"{name} − Analytical")
    ax_bot.set_xlabel(r"$\mathrm{time}~(\mu \mathrm{s})$")
    ax_bot.set_ylabel(r"$\Delta z~(\mu \mathrm{m})$")
    ax_bot.legend(loc="best")
    ax_bot.set_title(f"{name}: residuals")
    ax_bot.grid(True, linestyle="--", linewidth=0.8)

    fig.tight_layout()
    fig.savefig(outfile)

# ---- Produce figures ----
stacked_panels("RK4",   t,    z_rk, z_analytic_rk, "data/plot/rk4_z_and_residuals.pdf")
stacked_panels("Euler", t_eu, z_eu, z_analytic_eu, "data/plot/euler_z_and_residuals.pdf")

# plt.show()
plt.close("all")