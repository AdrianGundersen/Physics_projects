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
data1 = np.loadtxt("data/pos_vel_1.txt")
data2 = np.loadtxt("data/pos_vel_2.txt")

t1, x1, y1, z1, vx1, vy1, vz1 = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4], data1[:,5], data1[:,6]
t2, x2, y2, z2, vx2, vy2, vz2 = data2[:,0], data2[:,1], data2[:,2], data2[:,3], data2[:,4], data2[:,5], data2[:,6]

def plot_component(t1, arr1, t2, arr2, comp_name, xlabel, ylabel, filename):
    plt.figure()
    plt.plot(t1, arr1, alpha=0.8, label=f"Particle 1 {comp_name}")
    plt.plot(t2, arr2, alpha=0.8, label=f"Particle 2 {comp_name}")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


plot_component(t1, x1, t2, x2, "$x(t)$",
    r"$\text{time}~(\mu\text{s})$", r"$x\text{-position}~(\mu\text{m})$",
    "data/plot/x_position_vs_time.pdf")

plot_component(t1, y1, t2, y2, "$y(t)$",
    r"$\text{time}~(\mu\text{s})$", r"$y\text{-position}~(\mu\text{m})$",
    "data/plot/y_position_vs_time.pdf")

plot_component(t1, z1, t2, z2, "$z(t)$",
    r"$\text{time}~(\mu\text{s})$", r"$z\text{-position}~(\mu\text{m})$",
    "data/plot/z_position_vs_time.pdf")

plot_component(x1, vx1, x2, vx2, "$v_x(x)$",
    r"$x\text{-position}~(\mu\text{m})$", r"$v_x~(\mu\text{m}/\mu\text{s})$",
    "data/plot/v_x_vs_x.pdf")

plot_component(z1, vz1, z2, vz2, "$v_z(z)$",
    r"$z\text{-position}~(\mu\text{m})$", r"$v_z~(\mu\text{m}/\mu\text{s})$",
    "data/plot/v_z_vs_z.pdf")
