import numpy as np
import matplotlib.pyplot as plt
import os

# filenames: data/trapped_w0.200000-2.500000_dw0.005000_N100000.txt
w_min = 0.200000
w_max = 2.500000
w_step = 0.005000
N = 40000
C = 0 # coulomb on/off 0/1


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


filepath1 = f"data/trapped_w{w_min:.6f}-{w_max:.6f}_dw{w_step:.6f}_N{N}_C{C}.txt" # double precision


filepath2 = f"data/trapped_w{w_min:.6f}-{w_max:.6f}_dw{w_step:.6f}_N{N}_C{C-1}.txt" # comment out if only one file
plot_both = True  # set False to plot only filepath1

def load_file(path):
    data = np.loadtxt(path, comments="#")
    with open(path) as f:
        header = f.readline().strip().split()
    labels = [h for h in header if "frac_f" in h]
    omega = data[:, 0]
    frac = data[:, 1:]
    return omega, frac, labels


# ----- Physical parameters (consistent with your trap setup) -----
q = 1.0
m = 40.078 # Ca+
z0 = 20.0
d = 500.0
V0 = 25.0e-3 * 9.64852558e7
B = 9.64852558e1
omega_z = np.sqrt(2.0 * q * V0 / (m * d**2))
omega_0 = q * B / m
disc = omega_0**2 - 2.0 * omega_z**2
sqrt_disc = np.sqrt(disc)
omega_plus = 0.5 * (omega_0 + sqrt_disc)
omega_minus = 0.5 * (omega_0 - sqrt_disc)

DEN_TOL = 1e-10  # avoid division by small analytic values in relative errors

sources = []
if os.path.exists(filepath1):
    sources.append((filepath1, ""))  # suffix for legend
if plot_both and os.path.exists(filepath2) and filepath2 != filepath1:
    sources.append((filepath2, " (alt)"))

markers = ["o", "s", "^", "D", "x", "v", "P", "*"]
linestyles = ["-", "--"]

if not sources:
    raise FileNotFoundError("No input files found. Check filepath1/filepath2.")

base_labels = None
for si, (path, suff) in enumerate(sources):
    omega, frac_data, labels = load_file(path)
    if base_labels is None:
        base_labels = labels
    for i, label in enumerate(labels):
        f_value = label.split("f")[-1]
        m = markers[i % len(markers)]
        ls = linestyles[si % len(linestyles)]
        plt.plot(omega, frac_data[:, i], linestyle=ls, marker=".", label=rf"$f={float(f_value):.1f}$" + suff)

plt.xlabel(r"$\omega_V$ [MHz]")
plt.ylabel("Fraction of trapped particles")
plt.grid(True, which="both", alpha=0.3, linewidth=0.6, linestyle="--")
plt.legend(ncol=1)
plt.tight_layout()

outname = "data/plot/full_trapped.pdf" if len(sources) == 1 else "data/plot/full_trapped_combined.pdf"
os.makedirs(os.path.dirname(outname), exist_ok=True)
plt.savefig(outname)
