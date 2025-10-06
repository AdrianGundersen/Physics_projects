import numpy as np
import matplotlib.pyplot as plt

filepath1 = "data/fraction_trapped_vs_omega_zoomed_at_0.5_coloumb_off.txt"
filepath2 = "data/fraction_trapped_vs_omega_zoomed_at_0.5_coloumb_on.txt"

data = np.loadtxt(filepath1, comments="#")
with open(filepath1) as f:
    header = f.readline().strip().split()

# keep only fraction columns
labels = [h for h in header if "frac_f" in h]

omega = data[:, 0]
frac_data = data[:, 1:]

plt.rcParams.update({
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
plt.title("Trapping fraction vs. drive frequency")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("data/plot/zoomed_at_1.4_coloumb_off.pdf")
plt.show()
