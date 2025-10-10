import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 15,
    'figure.figsize': (6, 4),
    'axes.titlesize': 17,
    'axes.labelsize': 15,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'lines.linewidth': 2.0,
    'legend.fontsize': 11,
    'figure.dpi': 300,
})

os.makedirs("data/plot", exist_ok=True)

files = sorted(
    glob.glob("data/single_particle_N*.txt"),
    key=lambda p: int(re.search(r"_N(\d+)\.txt$", p).group(1))
)
if not files and os.path.exists("data/single_particle.txt"):
    files = ["data/single_particle.txt"]

# ----- Physical parameters (consistent with your trap setup) -----
q = 1.0
m = 1.0
z0 = 20.0
d = 500.0
V0 = 25.0e-3 * 9.64852558e7
B = 9.64852558e1
omega_z = np.sqrt(2.0 * q * V0 / (m * d**2))
omega_0 = q * B / m
DEN_TOL = 1e-10  # avoid division by small analytic values in relative errors

# ----- Analytic xy-solution builder; returns (x(t), y(t)) -----
def f(t):
    # initial conditions for radial motion
    x0 = 20.0
    v_0y = 25.0

    disc = omega_0**2 - 2.0 * omega_z**2
    sqrt_disc = np.sqrt(disc)
    omega_plus = 0.5 * (omega_0 + sqrt_disc)
    omega_minus = 0.5 * (omega_0 - sqrt_disc)

    # complex-amplitude solution for f(t) = x + i y
    A_plus  = (v_0y + omega_minus * x0) / (omega_minus - omega_plus)
    A_minus = -(v_0y + omega_plus  * x0) / (omega_minus - omega_plus)

    f_t = A_plus * np.exp(-1j * (omega_plus  * t)) + \
          A_minus* np.exp(-1j * (omega_minus * t))

    x = f_t.real
    y = f_t.imag
    return x, y

def short_N(n):
    # For naming eg. 32k
    n = int(n)
    return f"{n//1000}k" if n % 1000 == 0 else str(n) # if n < 1000, returns n

# ---- Figures ----
fig_rk, ax_rk = plt.subplots()
fig_eu, ax_eu = plt.subplots()
fig_rk_rel_r, ax_rk_rel_r = plt.subplots()
fig_eu_rel_r, ax_eu_rel_r = plt.subplots() 
fig_xy_rk, ax_xy_rk = plt.subplots() # plotting xy-plane for rk


for path in files:
    mN = re.search(r"_N(\d+)\.txt$", path) # searches for files matching the naming convention
    Nsuf = short_N(mN.group(1)) if mN else "N" # suffix eg. 32k

    data = np.loadtxt(path)
    # Columns: [0 t_rk, 1 x_rk, 2 y_rk, 3 z_rk, 4 vx_rk, 5 vy_rk, 6 vz_rk,
    #           7 t_eu, 8 x_eu, 9 y_eu, 10 z_eu, 11 vx_eu, 12 vy_eu, 13 vz_eu]

    # RK-data
    t = data[:, 0] # time-array
    x_rk, y_rk, z_rk = data[:, 1], data[:, 2], data[:, 3] # position of rk
    r_rk_mat = np.vstack([x_rk, y_rk, z_rk]) # shape = (3, N)
    r_rk = np.linalg.norm(r_rk_mat, axis=0) # shape = (N)


    # Forward Euler Data
    t_eu = data[:, 7]
    x_eu, y_eu, z_eu = data[:, 8], data[:, 9], data[:, 10]
    r_eu_mat = np.vstack([x_eu, y_eu, z_eu])
    r_eu = np.linalg.norm(r_eu_mat, axis=0)

    # analytical solutions
    z_an_rk = z0 * np.cos(omega_z * t) # z for rk-timesteps
    z_an_eu = z0 * np.cos(omega_z * t_eu) # --||-- forward euler
    x_an_rk, y_an_rk = f(t) # xy-plane for RK
    x_an_eu, y_an_eu = f(t_eu) 

    r_an_rk_mat = np.vstack([x_an_rk, y_an_rk, z_an_rk])
    r_an_rk = np.linalg.norm(r_an_rk_mat, axis=0)

    r_an_eu_mat = np.vstack([x_an_eu, y_an_eu, z_an_eu])
    r_an_eu = np.linalg.norm(r_an_eu_mat, axis=0)


    ax_xy_rk.plot(x_rk, y_rk, alpha=0.8, label=f"RK4 {Nsuf}")
    
    ax_rk.plot(t, z_rk, alpha=0.9, label=f"RK4 {Nsuf}")
    ax_eu.plot(t_eu, z_eu, alpha=0.9, label=f"EU {Nsuf}")


    mask_rk_r = r_an_rk > DEN_TOL
    mask_eu_r = r_an_eu > DEN_TOL
    rk_rel_r = np.abs(r_rk[mask_rk_r] - r_an_rk[mask_rk_r]) / r_an_rk[mask_rk_r]
    eu_rel_r = np.abs(r_eu[mask_eu_r] - r_an_eu[mask_eu_r]) / r_an_eu[mask_eu_r]

    ax_rk_rel_r.plot(t[mask_rk_r], rk_rel_r, alpha=0.5, label=f"RK4 {Nsuf}")
    ax_eu_rel_r.plot(t_eu[mask_eu_r], eu_rel_r, alpha=0.5, label=f"EU {Nsuf}")

    print(f"max relative error forward euler {path} = {np.max(eu_rel_r)} ")
    print(f"max relative error RK4           {path} = {np.max(rk_rel_r)} ")

# Add analytic z(t) reference to the z-plots for the densest N
if files:
    data_ref = np.loadtxt(files[-1])
    t_ref, t_eu_ref = data_ref[:, 0], data_ref[:, 7]
    ax_rk.plot(t_ref, z0 * np.cos(omega_z * t_ref), linestyle="--", alpha=0.7, label="Analytic")
    ax_eu.plot(t_eu_ref, z0 * np.cos(omega_z * t_eu_ref), linestyle="--", alpha=0.7, label="Analytic")
    x_an_ref, y_an_ref = f(t_ref)
    ax_xy_rk.plot(x_an_ref, y_an_ref, linestyle="--", alpha=0.8, color="black", label="Analytic")

# ---- Decorate & save ----
"""
ax_rk.set_title("RK4: z(t)")
ax_rk.set_xlabel(r"$t~(\mu \mathrm{s})$")
ax_rk.set_ylabel(r"$z~(\mu \mathrm{m})$")
ax_rk.legend(ncol=3)
fig_rk.tight_layout()
fig_rk.savefig("data/plot/rk4_z_vs_time_multiN.pdf")

ax_eu.set_title("Euler: z(t)")
ax_eu.set_xlabel(r"$t~(\mu \mathrm{s})$")
ax_eu.set_ylabel(r"$z~(\mu \mathrm{m})$")
ax_eu.legend(ncol=3)
fig_eu.tight_layout()
fig_eu.savefig("data/plot/euler_z_vs_time_multiN.pdf")
"""
ax_rk_rel_r.set_title(r"RK4: relative error $|r-r_a|/|r_a|$")
ax_rk_rel_r.set_xlabel(r"$t~(\mu \mathrm{s})$")
ax_rk_rel_r.set_ylabel(r"$|r-r_a|/|r_a|$")
ax_rk_rel_r.legend(ncol=3)
fig_rk_rel_r.tight_layout()
fig_rk_rel_r.savefig("data/plot/rk4_relerr_r_multiN.pdf")

ax_rk_rel_r.set_title(r"RK4: relative error $|r-r_a|/|r_a|$")
ax_rk_rel_r.set_xlabel(r"$t~(\mu \mathrm{s})$")
ax_rk_rel_r.set_ylabel(r"$\log\left(|r-r_a|/|r_a|\right)$")
ax_rk_rel_r.set_yscale("log")
ax_rk_rel_r.legend(ncol=3)
fig_rk_rel_r.tight_layout()
fig_rk_rel_r.savefig("data/plot/rk4_relerr_r_multiN_log.pdf")

ax_eu_rel_r.set_title(r"Euler: relative error $|r-r_a|/|r_a|$")
ax_eu_rel_r.set_xlabel(r"$t~(\mu \mathrm{s})$")
ax_eu_rel_r.set_ylabel(r"$|r-r_a|/|r_a|$")
ax_eu_rel_r.legend(ncol=3)
fig_eu_rel_r.tight_layout()
fig_eu_rel_r.savefig("data/plot/euler_relerr_r_multiN.pdf")

# RK4 x–y
ax_xy_rk.set_title(r"RK4: $x$–$y$ trajectory")
ax_xy_rk.set_xlabel(r"$x~(\mu \mathrm{m})$")
ax_xy_rk.set_ylabel(r"$y~(\mu \mathrm{m})$")
ax_xy_rk.set_xlim(-1,1)
ax_xy_rk.set_ylim(-19.8,-20.6)
ax_xy_rk.legend(ncol=2)
fig_xy_rk.tight_layout()
fig_xy_rk.savefig("data/plot/rk4_xy_multiN.pdf")


# plt.close("all")
# plt.show()