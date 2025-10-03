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
data = np.loadtxt("data/pos_vel.txt")

time = data[0,:]
x_position = data[1, :]
y_position = data[2, :]
z_position = data[3, :]

x_velocity = data[4, :]

plt.plot(time, y_position, label="z-position of particle 1")
plt.xlabel("time")
plt.ylabel("z-position")
plt.savefig("data/plot/z_position_vs_time.pdf")
