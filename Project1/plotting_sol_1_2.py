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
file_path = os.path.join(script_dir, "diff_eq_sol.txt")

data = np.loadtxt(file_path)

x = data[:, 0]
u = data[:, 1]

plt.plot(x, u)
plt.xlabel("x")
plt.ylabel("u(x)")
plt.xlim(min(x), max(x)*1.1)
plt.ylim(min(u), max(u)*1.1)
plt.title("Solution to the Differential Equation")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(file_path + "/Plot2")
plt.show()
plt.close()