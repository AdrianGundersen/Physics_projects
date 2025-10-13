import numpy as np
import matplotlib.pyplot as plt

# filenames: data/trapped_w0.200000-2.500000_dw0.005000_N100000.txt
w_min = 0.200000
w_max = 2.500000
w_step = 0.005000
N = 100000



filepath1 = f"data/trapped_w{w_min:.6f}-{w_max:.6f}_dw{w_step:.6f}_N{N}.txt" # double precision
print(filepath1)

"""
w_min = 0.1
w_max = 0.5
w_step = 0.05
N = 100
filepath2 = "data/fraction_trapped_vs_omega_zoomed_at_0.5_coloumb_on.txt"
"""

data = np.loadtxt(filepath1, comments="#")
with open(filepath1) as f:
    header = f.readline().strip().split()

# keep only fraction columns
labels = [h for h in header if "frac_f" in h]

omega = data[:, 0]
frac_data = data[:, 1:]

plt.rcParams.update({
    'font.family': 'serif',
    "font.size": 14,
    "figure.figsize": (6, 4),
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 2,
    "legend.fontsize": 10,
    "figure.dpi": 300,
})

for i, label in enumerate(labels):
    f_value = label.split("f")[-1]
    plt.plot(omega, frac_data[:, i], marker=".", label=rf"$f = {float(f_value):.1f}$")

plt.xlabel(r"$\omega_V$ [MHz]")
plt.ylabel("Fraction of trapped particles")
# plt.title("Trapping fraction vs. drive frequency")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("data/plot/zoomed_at_1.4_coloumb_off.pdf")
plt.show()
