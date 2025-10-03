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

t1, x1, y1, z1 = data1[:,0], data1[:,1], data1[:,2], data1[:,3]
t2, x2, y2, z2 = data2[:,0], data2[:,1], data2[:,2], data2[:,3]

def plot_component(t1, arr1, t2, arr2, comp_name, ylabel, filename):
    plt.figure()
    plt.plot(t1, arr1, label=f"Particle 1 {comp_name}(t)")
    plt.plot(t2, arr2, label=f"Particle 2 {comp_name}(t)")
    plt.xlabel("time (s)")
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

plot_component(t1, x1, t2, x2, "x", "x-position", "data/plot/x_position_vs_time.pdf")
plot_component(t1, y1, t2, y2, "y", "y-position", "data/plot/y_position_vs_time.pdf")
plot_component(t1, z1, t2, z2, "z", "z-position", "data/plot/z_position_vs_time.pdf")
