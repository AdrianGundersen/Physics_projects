# plotting_sol_1_2.py
# General plotting script for x against u(x) from a data file

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

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir)

folder = "output"
data_file_1 = "problem_8_1000.txt"
data_file_2 = "problem_8_100000.txt"

data_1 = np.loadtxt(os.path.join(folder, data_file_1))
data_2 = np.loadtxt(os.path.join(folder, data_file_2))

x_1, x_2 = data_1[:, 0], data_2[:, 0]
v_1, v_2 = data_1[:, 1], data_2[:, 1]
u_1, u_2 = data_1[:, 2], data_2[:, 2]
delta_1, delta_2 = data_1[:, 3], data_2[:, 3]
epsilon_1, epsilon_2 = data_1[:, 4], data_2[:, 4]

x_1, x_2 = x_1[1:-1], x_2[1:-1]
v_1, v_2 = v_1[1:-1], v_2[1:-1]
u_1, u_2 = u_1[1:-1], u_2[1:-1]
delta_1, delta_2 = delta_1[1:-1], delta_2[1:-1]
epsilon_1, epsilon_2 = epsilon_1[1:-1], epsilon_2[1:-1]

title_absolute = "Logaritmical absolute error"
title_relative = "Logaritmical relative error"

plt.plot(x_1, delta_1, label="n=1000")
plt.plot(x_2, delta_2, label="n=100000")
plt.xlabel("x")
plt.ylabel("log_10(|u-v|)")
plt.title(title_absolute)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(file_path + title_absolute + ".pdf")
plt.show()
plt.close()

plt.plot(x_1, epsilon_1, label="n=1000")
plt.plot(x_2, epsilon_2, label="n=100000")
plt.xlabel("x")
plt.ylabel("log_10(|(u-v)/u|)")
plt.title(title_relative)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(file_path + title_relative + ".pdf")
plt.show()
plt.close()