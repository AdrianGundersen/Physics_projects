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
    "1": 9.20388,
    "2": 4.64450,
    "4": 2.42251,
    "8": 1.64002,
    "16": 1.20351,
}

# cores and speedup
x = np.array([int(k) for k in runtimes.keys()], dtype=float)  # number of cores
y = np.array([v for v in runtimes.values()], dtype=float)     # runtimes
speedup = y[0] / y                                            # S(p) = T1 / Tp


# ideal scaling S_ideal(p) = p
ideal_speedup = x

plt.figure()
plt.plot(x, speedup, 'o', label='Measured speedup')
plt.plot(x, ideal_speedup, ':', label='Ideal speedup')
plt.grid()
plt.xlabel('Number of cores')
plt.ylabel('Speedup $S = T_1 / T_p$')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "speedup_plot.png", bbox_inches='tight')
plt.close()