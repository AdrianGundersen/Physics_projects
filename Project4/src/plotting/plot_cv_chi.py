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

fig_dir = ROOT / "data/figures/Cv_chi"
fig_dir.mkdir(parents=True, exist_ok=True)


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
    plt.savefig(fig_dir / f'm_vs_T_L{L}.pdf')
    plt.close()
    return None

def Tc_regress(data, L, observable="Cv", plot = False, min_sweeps=1e4):
    """
    """
    T_values = []
    obs_values = []

    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        if values["sweeps"] < min_sweeps:
            continue # skip if not enough sweeps
        if observable == "Cv":
            obs_values.append(values["Cv"])
        elif observable == "chi":
            obs_values.append(values["chi"])
        else:
            raise ValueError("Observable must be 'Cv' or 'chi'.")

    obs_max = max(obs_values)
    max_index = obs_values.index(obs_max)
    delta_idx = 5 # number of points on each side of the maximum to consider
    start_idx = max(0, max_index - delta_idx)
    end_idx = min(len(T_values), max_index + delta_idx + 1)
    T_fit = np.array(T_values[start_idx:end_idx])
    obs_fit = np.array(obs_values[start_idx:end_idx])
    coeffs = np.polyfit(T_fit, obs_fit, 2)
    poly = np.poly1d(coeffs)
    T_c = -coeffs[1] / (2 * coeffs[0]) # -b/2a
    observable_max = poly(T_c) 

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
        plt.grid()
        plt.savefig(fig_dir / f'{observable}_Tc_fit_L{L}.pdf')
        plt.close()

    return T_c, observable_max
# Load data and plot for L
L = np.array([10, 20])#, 60, 70, 80, 90, 100, 110, 120, 130])
min_sweeps = 1e5
json_path = ROOT / "test.json"
raw_data = load_JSON(json_path)
data = sort_data(raw_data, min_sweeps=min_sweeps)

Tc_Cv_vals = []
Tc_chi_vals = []

for l in L:
    plot_Cv_vs_T(data, l)
    plot_chi_vs_T(data, l)
    plot_eps_vs_T(data, l)
    plot_m_vs_T(data, l)
    Tc_Cv, Cv_max = Tc_regress(data, l, observable="Cv", plot=True, min_sweeps=min_sweeps)
    Tc_chi, chi_max = Tc_regress(data, l, observable="chi", plot=True, min_sweeps=min_sweeps)
    Tc_Cv_vals.append(Tc_Cv)    
    Tc_chi_vals.append(Tc_chi)
    print(f"L={l}: Tc from Cv = {Tc_Cv:.4f} (max Cv = {Cv_max:.4f}), Tc from chi = {Tc_chi:.4f} (max chi = {chi_max:.4f})")

# 1/L vs Tc plot
Tc_Cv_vals = np.array(Tc_Cv_vals)
Tc_chi_vals = np.array(Tc_chi_vals)
Tc_vals = 0.5 * (Tc_Cv_vals + Tc_chi_vals)

inv_L = 1.0 / L

# linregress
slope, intercept, r_value, p_val, std_err = linregress(inv_L, Tc_vals)


plt.figure()
plt.plot(inv_L, Tc_vals, 'o', label='Tc')
plt.plot(inv_L, intercept + slope * inv_L, 'r--', label=f'Fit Cv: Tc={intercept:.3f} + {slope:.3f}/L')
plt.xlabel('1/L')
plt.ylabel('Critical Temperature Tc')
plt.legend()
plt.grid()
plt.savefig(fig_dir / 'Tc_vs_invL.pdf')
plt.close()