# Heavily chatGPT made plots as mentioned in report.
import os
import numpy as np
import matplotlib.pyplot as plt

# ---- Matplotlib defaults ----
plt.rcParams.update({
    'font.size': 14,
    'figure.figsize': (12, 4),  # two panels side-by-side
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2,
    'legend.fontsize': 10,
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

def double_panel(name: str, t_num, z_num, z_an, outfile: str):
    """
    Left: numerical vs analytical z(t).
    Right: residuals (numerical - analytical).
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # Left: trajectories
    ax = axes[0]
    ax.plot(t_num, z_num, label=f"{name} z(t)")
    ax.plot(t_num, z_an, "--", label="Analytical z(t)")
    ax.set_xlabel(r"$\mathrm{time}~(\mu s)$")
    ax.set_ylabel(r"$z~(\mu m)$")
    ax.legend()
    ax.set_title(f"{name}: z(t) vs analytical")

    # Right: residuals
    ax = axes[1]
    res = z_num - z_an
    ax.plot(t_num, res, label=f"{name} − Analytical")
    ax.set_xlabel(r"$\mathrm{time}~(\mu s)$")
    ax.set_ylabel(r"$\Delta z~(\mu m)$")
    ax.legend()
    ax.set_title(f"{name}: residuals")

    fig.tight_layout()
    fig.savefig(outfile)

# ---- Produce figures ----
double_panel("RK4",  t,    z_rk, z_analytic_rk, "data/plot/rk4_z_and_residuals.pdf")
double_panel("Euler", t_eu, z_eu, z_analytic_eu, "data/plot/euler_z_and_residuals.pdf")

plt.show()
