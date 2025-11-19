# plotting/plot_speedup.py
import numpy as np
import matplotlib as mpl
mpl.use("Agg")  # to avoid wayland issues
import matplotlib.pyplot as plt
from scipy.stats import linregress
from pathlib import Path

# plotting style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 14,
    'figure.figsize': (6, 4),
    'axes.titlesize': 14,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.0,
    'legend.fontsize': 14,
    'figure.dpi': 300,
})

# Finds project root
ROOT = Path(__file__).resolve().parents[2]

fig_dir = ROOT / "Project4/data/figures/"
fig_dir.mkdir(parents=True, exist_ok=True)

runtimes = {  # manually inputted
    "1": 25.9618,
    "2": 13.1849,
    "3": 8.94716,
    "4": 6.81619,
    "6": 4.55323,
    "8": 3.51539,
    "12": 2.46687,
    "16": 2.02883,
}

# cores and speedup
x = np.array([int(k) for k in runtimes.keys()], dtype=float)  # number of cores
y = np.array([v for v in runtimes.values()], dtype=float)     # runtimes
speedup = y[0] / y                                            # S(p) = T1 / Tp


# ideal scaling S_ideal(p) = p
ideal_speedup = x

time_saved_factor = runtimes["1"]/runtimes["12"]
print(f"Factor for saved trime = {time_saved_factor:.1f}")

plt.figure()
plt.plot(x, ideal_speedup, '--', label='Ideal speedup')
plt.plot(x, speedup, 'o', label='Measured speedup')
plt.grid(alpha=0.6, linestyle='--')
plt.xlabel('Number of threads')
plt.ylabel('Speedup $t_1 / t_p$')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "speedup_plot.pdf", bbox_inches='tight')
plt.close()
