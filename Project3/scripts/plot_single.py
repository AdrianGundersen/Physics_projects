# Heavily chatGPT made plots as mentioned in report.
import os
import numpy as np
import matplotlib.pyplot as plt

# Oppdaterer mpl parametere
plt.rcParams.update({
    'font.size': 14,
    'figure.figsize': (6, 4),
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2,
    'legend.fontsize': 10,
    'figure.dpi': 300,
})

os.makedirs("data/plot", exist_ok=True)
data1 = np.loadtxt("data/single_particle.txt")

t1, x1, y1, z1, vx1, vy1, vz1 = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4], data1[:,5], data1[:,6]
V0_over_d2 = 9.65
q = 1
ca_mass = 40.078
z_0 = 20
d  = 500.0
V0 = 2.41e6
z1 = z1

# multiply by 2pi to get rad
omega_z = np.sqrt(2*q*V0/(ca_mass*d**2)) 
phi = omega_z * (2*np.pi) 
z_t = z_0*np.cos(omega_z*t1)

def plot_component(t1, arr1, t2, arr2, comp_name, xlabel, ylabel, filename):
    plt.figure()
    plt.plot(t1, arr1, alpha=0.8, label=f"Simulated {comp_name}")
    plt.plot(t2, arr2, alpha=0.8, label=f"Analytical {comp_name}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

plot_component(t1, z1, t1, z_t, "$z(t)$",
    r"$\text{time}~(\mu\text{s})$", r"$z\text{-position}~(\mu\text{m})$",
    "data/plot/z_t_vs_analytical.pdf")
