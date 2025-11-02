import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# plotting style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 15,
    'figure.figsize': (6, 4),
    'axes.titlesize': 17,
    'axes.labelsize': 20,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.0,
    'legend.fontsize': 17,
    'figure.dpi': 300,
})

# Finds project root
ROOT = Path(__file__).resolve().parents[2]

fig_dir = ROOT / "data/figures/L20"
fig_dir.mkdir(parents=True, exist_ok=True)

def load_data(T, spin_config):
    p = ROOT / f"data/outputs/L20_T={T:0.6f}_spin={spin_config}.txt"
    if not p.exists():
        print(f"Warning: missing file: {p}")
        return None
    try:
        return np.loadtxt(p)
    except Exception as e:
        print(f"Error reading {p}: {e}")
        return None
    
def cummean(x):
    return np.cumsum(x) / np.arange(1, len(x) + 1)


# Spin configs and temps
spin_configs = ["all_up", "all_down", "random"]
temperatures = [1.0, 2.40]  # flere kan legges til

epsT10_u = load_data(1.0, "random")
epsT10_o = load_data(1.0, "all_up")
epsT24_u = load_data(2.40, "random")
epsT24_o = load_data(2.40, "all_up")
mean_eps_T10_u = cummean(epsT10_u)
mean_eps_T10_o = cummean(epsT10_o)
mean_eps_T24_u = cummean(epsT24_u)
mean_eps_T24_o = cummean(epsT24_o)

n_MC_cycles = np.arange(0, len(epsT10_u))

plt.plot(n_MC_cycles, epsT10_u, '-', color='#377eb8', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, epsT10_o, '-', color='#4daf4a', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, epsT24_u, '-', color='#e41a1c', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, epsT24_o, '-', color='#984ea3', alpha=0.4, linewidth=1.0)
plt.plot(n_MC_cycles, mean_eps_T10_u, '-', linewidth=2.0, color='#377eb8', label='$T=1.0$ $J/k_{B}$, unordered')
plt.plot(n_MC_cycles, mean_eps_T10_o, '-', linewidth=2.0, color='#4daf4a', label='$T=1.0$ $J/k_{B}$, ordered')
plt.plot(n_MC_cycles, mean_eps_T24_u, '-', linewidth=2.0, color='#e41a1c', label='$T=2.4$ $J/k_{B}$, unordered')
plt.plot(n_MC_cycles, mean_eps_T24_o, '-', linewidth=2.0, color='#984ea3', label='$T=2.4$ $J/k_{B}$, ordered')
plt.xscale('log')
plt.xlabel("Monte Carlo cycles")
plt.ylabel(r"$\epsilon$")
plt.grid(True, which='both', alpha=0.25)
plt.legend(loc='best', frameon=False)

out = fig_dir / "burnin_L20_eps_and_running_mean.pdf"
plt.tight_layout()
plt.savefig(out, format="pdf", bbox_inches="tight")
plt.close()
print(f"Saved: {out}")