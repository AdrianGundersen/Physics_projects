# src/plotting/plot_Cv.py
import numpy as np
import matplotlib as mpl
mpl.use("Agg") # to avoid wayland issues
import matplotlib.pyplot as plt
import json
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

fig_dir = ROOT / "Project4/data/figures/Cv_chi"
fig_dir.mkdir(parents=True, exist_ok=True)

def get_max(data, L, observable="Cv"):
    """
    Get maximum value of observable and index for given L.
    """
    L_data = data[str(L)]
    obs_vals = [L_data[T_str][observable] for T_str in L_data]
    return max(obs_vals), obs_vals.index(max(obs_vals))

def uncertainty_cv_chi(data, L, index_cv, index_chi):
    """
    Calculate uncertainty in chi and Cv using variance of energy and magnetization.
    """

    L_data = data[str(L)]

    t_sorted = sorted(L_data.keys(), key=float)

    T_cv_str = t_sorted[index_cv]
    vals_cv = L_data[T_cv_str]
    T_cv = float(T_cv_str)

    T_chi_str = t_sorted[index_chi]
    vals_chi = L_data[T_chi_str]
    T_chi = float(T_chi_str)

    N = L * L
    kB = 1.0

    M_cv = vals_cv["sweeps"] * vals_cv["walkers"]
    M_chi = vals_chi["sweeps"] * vals_chi["walkers"]

    eps = vals_cv["avg_eps"]
    eps2 = vals_cv["avg_eps2"]
    eps4 = vals_cv["avg_eps4"]

    m = vals_chi["avg_mabs"]
    m2 = vals_chi["avg_mabs2"]
    m4 = vals_chi["avg_mabs4"]

    var_eps = eps2 - eps**2
    var_m = m2 - m**2
    var_eps2 = eps4 - eps2**2
    var_m2 = m4 - m2**2

    sigma_cv = N/(kB*T_cv**2)*np.sqrt(var_eps2/M_cv + 4*eps**2*var_eps/M_cv)
    sigma_chi = N/(kB*T_chi)*np.sqrt(var_m2/M_chi + 4*m**2*var_m/M_chi)

    return sigma_chi, sigma_cv

def load_JSON(filepath) -> dict:
    """Load JSON data from a file.
    
    structure of json: { L : { T : {chi, Cv, sweeps, walkers} } }
    """
    with open(filepath, 'r') as f:
        data = json.load(f) # loads into dictionary
    return data

def sort_data(data, min_sweeps=1e4):
    """
    Sort data dictionary by temperature for each L. Remove if sweeps < min_sweeps.
    """

    sorted_data = {}
    for L_str, T_dict in data.items():
        T_values = []
        for T_str, values in T_dict.items():
            if values["sweeps"] >= min_sweeps:
                T_values.append((float(T_str), values))
            T_values.sort(key=lambda x: x[0])  # sort by temperature
        sorted_data[L_str] = {str(T): vals for T, vals in T_values}
    return sorted_data


def plot_Cv_vs_T(data, L):
    """
    Plot heat capacity Cv vs temperature T for a given L.
    """
    T_values = []
    Cv_values = []

    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        Cv_values.append(values["Cv"])

    plt.figure()
    plt.plot(T_values, Cv_values, marker='o', linestyle='-')

    plt.xlabel('Temperature T')
    plt.ylabel('Heat Capacity Cv')
    plt.tight_layout()
    plt.grid()
    plt.savefig(fig_dir / f'Cv_vs_T_L{L}.pdf')
    plt.close()
    return None

def plot_chi_vs_T(data, L):
    """
    Plot magnetic susceptibility chi vs temperature T for a given L.
    """
    T_values = []
    chi_values = []

    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        chi_values.append(values["chi"])

    plt.figure()
    plt.plot(T_values, chi_values, marker='o', linestyle='-')
    plt.xlabel('Temperature T')
    plt.ylabel('Magnetic Susceptibility χ')
    plt.grid()
    plt.tight_layout()
    plt.savefig(fig_dir / f'chi_vs_T_L{L}.pdf')
    plt.close()
    return None

def plot_eps_vs_T(data, L):
    """
    Plot energy per spin ε vs temperature T for a given L.
    """
    T_values = []
    eps_values = []

    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        eps_values.append(values["avg_eps"])  # energy per spin

    plt.figure()
    plt.plot(T_values, eps_values, marker='o', linestyle='-')
    plt.xlabel('Temperature T')
    plt.ylabel('Energy per spin ε')
    plt.grid()
    plt.tight_layout()
    plt.savefig(fig_dir / f'eps_vs_T_L{L}.pdf')
    plt.close()
    return None

def plot_m_vs_T(data, L):
    """
    Plot absolute magnetization per spin m vs temperature T for a given L.
    """
    T_values = []
    m_values = []

    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        m_values.append(values["avg_mabs"])  # magnetization per spin

    plt.figure()
    plt.plot(T_values, m_values, marker='o', linestyle='-')
    plt.xlabel('Temperature T')
    plt.ylabel('Magnetization per spin m')
    plt.grid()
    plt.tight_layout()
    plt.savefig(fig_dir / f'm_vs_T_L{L}.pdf')
    plt.close()
    return None

def Tc_regress(data, L, observable="Cv", plot = False):
    """
    """
    L_data = data[str(L)]
    T_sorted = sorted(L_data.keys(), key=float)
    T_values = []
    obs_values = []

    for T_str in T_sorted:
        values = L_data[T_str]
        T_values.append(float(T_str))
        if observable == "Cv":
            obs_values.append(values["Cv"])
        elif observable == "chi":
            obs_values.append(values["chi"])
        else:
            raise ValueError("Observable must be 'Cv' or 'chi'.")

    obs_max = max(obs_values)
    max_index = obs_values.index(obs_max)
    delta_idx = 4 # number of points on each side of the maximum to consider
    start_idx = max(0, max_index - delta_idx)
    end_idx = min(len(T_values), max_index + delta_idx + 1)

    T_fit = np.array(T_values[start_idx:end_idx])
    obs_fit = np.array(obs_values[start_idx:end_idx])

    coeffs, cov = np.polyfit(T_fit, obs_fit, 2, cov=True)
    a, b, c = coeffs
    poly = np.poly1d(coeffs)

    T_c = -coeffs[1] / (2 * coeffs[0]) # -b/2a
    observable_max = poly(T_c) 
    
    var_a = cov[0,0]
    var_b = cov[1,1]
    var_ab = cov[0,1]

    dTda = b / (2.0 * a**2)
    dTdb = -1.0 / (2.0 * a)
    sigma_Tc = np.sqrt(dTda**2 * var_a + dTdb**2 * var_b + 2 * dTda * dTdb * var_ab)

    if observable == "Cv":
        sigma_obs_max = uncertainty_cv_chi(data, L, max_index, max_index)[1]
    elif observable == "chi":
        sigma_obs_max = uncertainty_cv_chi(data, L, max_index, max_index)[0]
    else:
        raise ValueError("Observable must be 'Cv' or 'chi'.")

    # Plotting
    if plot:
        plt.figure()
        T_plot = np.linspace(min(T_fit) - abs(min(T_fit)*0.01), max(T_fit) + abs(max(T_fit)*0.01), 100)
        plt.plot(T_plot, poly(T_plot), label='Fitted curve')
        plt.scatter(T_fit, obs_fit, color='orange', label='Data points')
        plt.axvline(T_c, color='r', linestyle='--', label=f'T_c = {T_c:.2f}')
        plt.xlabel('Temperature T')
        plt.ylabel(f'{observable} (max: {observable_max:.2f} at T_c: {T_c:.2f})')
        plt.legend()
        plt.tight_layout()
        plt.grid()
        plt.savefig(fig_dir / f'{observable}_Tc_fit_L{L}.pdf')
        plt.close()

    return T_c, observable_max, sigma_obs_max, sigma_Tc


def plot_together(data, L, observable="Cv"):
    """
    Plot observable for all L in the same plot
    """
    T_values = []
    obs_values = []

    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        obs_values.append(values[observable])

    plt.plot(T_values, obs_values, linestyle='-', label=f'L={L}', alpha = 0.7)
    plt.legend()
    plt.grid()
    return None



# Load data and plot for L
L = np.array([5, 10, 20, 30])# 40, 50, 100])#, 60, 70, 80, 90, 100, 110, 120, 130])
min_sweeps = 1e5
json_path = ROOT / "Project4/data/outputs/tsweep_results.json"
raw_data = load_JSON(json_path)
data = sort_data(raw_data, min_sweeps=min_sweeps)

Tc_Cv_vals = []
Tc_chi_vals = []
sigma_TcCv_vals = []
sigma_TcChi_vals = []
sigma_Cv_vals = []
sigma_chi_vals = []

for l in L:
    plot_Cv_vs_T(data, l)
    plot_chi_vs_T(data, l)
    plot_eps_vs_T(data, l)
    plot_m_vs_T(data, l)
    Tc_Cv, Cv_max, sigma_Cv, sigma_TcCv = Tc_regress(data, l, observable="Cv", plot=True)
    Tc_chi, chi_max, sigma_chi, sigma_TcChi = Tc_regress(data, l, observable="chi", plot=True)
    Tc_Cv_vals.append(Tc_Cv)
    Tc_chi_vals.append(Tc_chi)
    sigma_TcCv_vals.append(sigma_TcCv)
    sigma_TcChi_vals.append(sigma_TcChi)
    sigma_Cv_vals.append(sigma_Cv)
    sigma_chi_vals.append(sigma_chi)
    print(f"L={l}: Tc from Cv = {Tc_Cv:.4f}+-{sigma_TcCv:.4f} (max Cv = {Cv_max:.4f}+-{sigma_Cv:.4f}), Tc from chi = {Tc_chi:.4f}+-{sigma_TcChi:.4f} (max chi = {chi_max:.4f}+-{sigma_chi:.4f})")

# 1/L vs Tc plot
Tc_Cv_vals = np.array(Tc_Cv_vals)
Tc_chi_vals = np.array(Tc_chi_vals)
sigma_TcCv_vals = np.array(sigma_TcCv_vals)
sigma_TcChi_vals = np.array(sigma_TcChi_vals)

Tc_vals = 0.5 * (Tc_Cv_vals + Tc_chi_vals)
Tc_errs = 0.5 * np.sqrt(sigma_TcCv_vals**2 + sigma_TcChi_vals**2)

inv_L = 1.0 / L

wheights = 1.0 / Tc_errs**2

# linregress
slope, intercept, r_value, p_val, std_err = linregress(inv_L, Tc_vals)

N = len(L)
inv_L_mean = np.mean(inv_L)
intercept_err = std_err * np.sqrt(np.sum(inv_L**2) / (N*np.sum((inv_L - inv_L_mean)**2)))

print(f"Tc(inf) = {intercept:.5f} ± {intercept_err:.5f}")

plt.figure()
plt.plot(inv_L, Tc_vals, 'o', label='Tc')
plt.plot(inv_L, intercept + slope * inv_L, 'r--', label=f'Fit Cv: Tc={intercept:.3f} + {slope:.3f}/L')
plt.xlabel('1/L')
plt.ylabel('Critical Temperature Tc')
plt.legend()
plt.tight_layout()
plt.grid()
plt.savefig(fig_dir / 'Tc_vs_invL.pdf')
plt.close()

# plot together
plt.figure()
for observable in ["Cv", "chi", "avg_eps", "avg_mabs"]:
    plt.clf()
    for l in L:
        plot_together(data, l, observable)
    plt.xlabel('Temperature T')
    plt.ylabel(f'{observable}')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(fig_dir / f'{observable}_vs_T_all_L.pdf')
    plt.close()