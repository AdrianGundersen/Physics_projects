# plotting_difference.py
# Plots the absolute and relative difference between numerical and analytical solutions from a data file.

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
file_path = os.path.join(script_dir, "output")

n_vals = np.array([1e1, 1e2, 1e3, 1e4, 1e5])

analytical = np.loadtxt(os.path.join(file_path, "diff_eq_sol.txt"))
plt.figure()
plt.title("Analytical Solution vs Numerical Solutions")
plt.xlabel(r"$x$")
plt.ylabel(r"$u(x)$")
plt.plot(analytical[:,0], analytical[:,1], label="Analytical Solution", color='black', alpha = 1)


for n in n_vals:
    data_file = f"problem_8_{int(n)}.txt"
    data = np.loadtxt(os.path.join(file_path, data_file))
    x_vals = np.array(data[:,0])
    u_numerical = np.array(data[:,1])

    # Adding boundary points
    x_vals = np.insert(x_vals, 0, 0)
    x_vals = np.append(x_vals, 1)
    u_numerical = np.insert(u_numerical, 0, 0)
    u_numerical = np.append(u_numerical, 0)

    plt.plot(x_vals, u_numerical, label=f"Numerical n={int(n)}", alpha = 0.6, linestyle='--')


plt.legend(loc = "upper right")
plt.tight_layout()
plt.grid(alpha=0.3)
plt.savefig(os.path.join(file_path, "analytical_vs_numerical.pdf"))
plt.show()

