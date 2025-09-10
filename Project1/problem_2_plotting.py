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
data_file = input("Enter the data file name: ")
title = input("Enter the plot title: ")

data = np.loadtxt(os.path.join(folder, data_file))

x = data[:, 0]
u = data[:, 1]

plt.plot(x, u)
plt.xlabel(r"$x$")
plt.ylabel(r"$u(x)$")
plt.xlim(min(x), max(x)*1.1)
plt.ylim(min(u), max(u)*1.1)
plt.title(title)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(file_path, title) + ".pdf")
plt.show()
plt.close()