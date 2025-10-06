# plot_trapped_vs_omega.py
# Heavily chatGPT made plots as mentioned in report.
import numpy as np
import matplotlib.pyplot as plt

# ---- Load data ----
data = np.loadtxt("data/fraction_trapped_vs_omega.txt", comments="#")
with open("data/fraction_trapped_vs_omega.txt") as f:
    header = f.readline().strip().split()

# Extract labels and data
omega = data[:, 0]
frac_data = data[:, 1:]
labels = header[1:]  # e.g. ["frac_f0.100000", "frac_f0.400000", ...]

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

for i, label in enumerate(labels):
    f_value = label.split("f")[-1]
    plt.plot(omega, frac_data[:, i], marker="o", label=f"$f = {float(f_value):.1f}$")

plt.xlabel(r"$\omega_V$ [MHz]")
plt.ylabel("Fraction of trapped particles")
plt.title("Trapping fraction vs. drive frequency")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
