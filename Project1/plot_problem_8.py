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

data = np.loadtxt(os.path.join(folder, data_file))

x = data[:, 0]
v = data[:, 1]
u = data[:, 2]
delta = data[:, 3]

x = x[1:-1]
v = v[1:-1]
u = u[1:-1]
delta = delta[1:-1]

title = "u(x) vs. v(x)"

plt.plot(x, v, label="v(x)")
plt.plot(x, u, label="u(x)", linestyle="-")
plt.xlabel("x")
plt.ylabel("u, v")
plt.title(title)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(file_path + title + ".pdf")
plt.show()
plt.close()

title_log = "Logaritmical absolute error"
plt.plot(x, delta, label="v(x)")
plt.xlabel("x")
plt.ylabel("log_10(|u-v|)")
plt.title(title_log)
plt.grid(alpha=0.3)
#plt.legend()
plt.tight_layout()
plt.savefig(file_path + title_log + ".pdf")
plt.show()
plt.close()