import numpy as np
import matplotlib as mpl
mpl.use("Agg") # to avoid wayland issues
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
    """_summary_

    Args:
        T (_type_): _description_
        spin_config (_type_): _description_

    Returns:
        _type_: _description_
    """    
    p = ROOT / f"data/outputs/L20_T={T:0.6f}_spin={spin_config}.txt"
    data = np.genfromtxt(
        p,
        delimiter=',',
        comments='#',      # skip "# m, e"
        usecols=(0,1),    # only load energy and magnetization
        autostrip=True,
        invalid_raise=False
    )
    if data.ndim == 1:     # handle single-row files
        data = data.reshape(1, -1)
    return data


def cummean(x):
    return np.cumsum(x) / np.arange(1, len(x) + 1)


# Spin configs and temps
spin_configs = ["all_up", "all_down", "random"]
temperatures = [1.0, 2.40]  # flere kan legges til

epsT10_u = load_data(1.0, "random")[:, 0]
epsT10_o = load_data(1.0, "all_up")[:, 0]
epsT24_u = load_data(2.40, "random")[:, 0]
epsT24_o = load_data(2.40, "all_up")[:, 0]
mean_eps_T10_u = cummean(epsT10_u)
mean_eps_T10_o = cummean(epsT10_o)
mean_eps_T24_u = cummean(epsT24_u)
mean_eps_T24_o = cummean(epsT24_o)

n10u = np.arange(len(epsT10_u))
n10o = np.arange(len(epsT10_o))
n24u = np.arange(len(epsT24_u))
n24o = np.arange(len(epsT24_o))

print(rf"length of different data sets: {len(epsT10_u)}, {len(epsT10_o)}, {len(epsT24_u)}, {len(epsT24_o)}")

plt.plot(n10u, epsT10_u, '-', color='#377eb8', alpha=0.4, linewidth=1.0)
plt.plot(n10o, epsT10_o, '-', color='#4daf4a', alpha=0.4, linewidth=1.0)
plt.plot(n24u, epsT24_u, '-', color='#e41a1c', alpha=0.4, linewidth=1.0)
plt.plot(n24o, epsT24_o, '-', color='#984ea3', alpha=0.4, linewidth=1.0)

plt.plot(n10u, mean_eps_T10_u, '-', linewidth=2.0, color='#377eb8', label='$T=1.0$ $J/k_{B}$, unordered')
plt.plot(n10o, mean_eps_T10_o, '-', linewidth=2.0, color='#4daf4a', label='$T=1.0$ $J/k_{B}$, ordered')
plt.plot(n24u, mean_eps_T24_u, '-', linewidth=2.0, color='#e41a1c', label='$T=2.4$ $J/k_{B}$, unordered')
plt.plot(n24o, mean_eps_T24_o, '-', linewidth=2.0, color='#984ea3', label='$T=2.4$ $J/k_{B}$, ordered')

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

# plot histograms after burn-in

burn_in = 1e4 # discard first 10 000 sweeps
burnt_in_epsT10_u = epsT10_u[int(burn_in):]
burnt_in_epsT24_u = epsT24_u[int(burn_in):]

def plot_hist_gauss(data, T, fig_dir):
    """Plot histogram and Gaussian fit for the given data.

    Args:
        data (np.ndarray): Data to plot.
        T (float): Temperature parameter.
        fig_dir (Path): Directory to save the figure.

    Returns:
        None
    """    
    plt.figure()

    bin_width = 0.01 # 4J / 400 particles delta epsilon
    plt.hist(data, bins=np.arange(min(data), max(data) + bin_width, bin_width), density=True, alpha=0.6, color='g', label='Histogram')

    mu = np.mean(data)
    sigma = np.std(data)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    plt.plot(x, p, 'k', linewidth=2, label='Gaussian fit')
    plt.axvline(mu, color='r', linestyle='dashed', linewidth=1, label=r'$\mu$')
    plt.axvline(mu + sigma, color='b', linestyle='dashed', linewidth=1, label=r'$\mu \pm \sigma$')
    plt.axvline(mu - sigma, color='b', linestyle='dashed', linewidth=1)
    plt.xlim(xmin, xmax)
    plt.xlabel(r'$\epsilon$')
    plt.ylabel('Density')
    plt.legend()
    
    out = fig_dir / f"burnin_L20_hist_eps_T={T:0.2f}.pdf"
    plt.tight_layout()
    plt.savefig(out, format="pdf", bbox_inches="tight")
    plt.close()
    print(rf"Mean: {mu:.4f}, Std: {sigma:.4f} for T={T:0.2f}")
    print(f"Saved: {out}")


plot_hist_gauss(burnt_in_epsT10_u, 1.0, fig_dir)
plot_hist_gauss(burnt_in_epsT24_u, 2.4, fig_dir)

