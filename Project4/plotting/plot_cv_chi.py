# src/plotting/plot_Cv.py
# With much help from ChatGPT
import numpy as np
import matplotlib as mpl
mpl.use("Agg") # to avoid wayland issues
import matplotlib.pyplot as plt
import json
from scipy.stats import linregress
import math
from pathlib import Path

# plotting style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 14,
    'figure.figsize': (6, 4),
    'axes.titlesize': 14,
    'axes.labelsize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.0,
    'legend.fontsize': 15,
    'figure.dpi': 300,
})

# Finds project root
ROOT = Path(__file__).resolve().parents[2]

fig_dir = ROOT / "Project4/data/figures/Cv_chi"
fig_dir.mkdir(parents=True, exist_ok=True)
L = np.array([5, 7, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120])
min_sweeps = 1e5 + 1
min_sweeps_fit = 1e5
json_path = ROOT / "Project4/data/outputs/tsweep_finale.json"

plot_all = False # whether to plot

def load_JSON(filepath) -> dict:
    """Load JSON data from a file.
    
    structure of json: { L : { T : {chi, Cv, sweeps, walkers} } }
    """
    with open(filepath, 'r') as f:
        data = json.load(f) # loads into dictionary
    return data


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
    eps3 = vals_cv["avg_eps3"]
    eps4 = vals_cv["avg_eps4"]

    m = vals_chi["avg_mabs"]
    m2 = vals_chi["avg_mabs2"]
    m3 = vals_chi["avg_mabs3"]
    m4 = vals_chi["avg_mabs4"]

    var_eps = eps2 - eps**2
    var_m = m2 - m**2
    var_eps2 = eps4 - eps2**2
    var_m2 = m4 - m2**2

    cov_eps = eps3 - eps*eps2
    cov_m = m3 - m*m2

    sigma_cv = N/(kB*T_cv**2)*np.sqrt(var_eps2/M_cv + 4*eps**2*var_eps/M_cv - 4*eps* cov_eps/M_cv)
    sigma_chi = N/(kB*T_chi)*np.sqrt(var_m2/M_chi + 4*m**2*var_m/M_chi - 4*m*cov_m/M_chi)

    return sigma_chi, sigma_cv

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
    plt.ylabel(r'Heat Capacity $c_v$')
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
    plt.ylabel(r'Magnetic Susceptibility $\chi_s$')
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
    plt.ylabel(r'Energy per spin $\langle \epsilon \rangle$')
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
    plt.ylabel(r'Magnetization per spin $\langle|m|\rangle$')
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
    delta_idx = 3 # number of points on each side of the maximum to consider
    start_idx = max(0, max_index - delta_idx)
    end_idx = min(len(T_values), max_index + delta_idx + 1)

    sigma_obs_vals = []
    for i, T_str in enumerate(T_sorted):
        if observable == "Cv":
            sigma_obs = uncertainty_cv_chi(data, L, i, i)[1]
        elif observable == "chi":
            sigma_obs = uncertainty_cv_chi(data, L, i, i)[0]
        else:
            raise ValueError("Observable must be 'Cv' or 'chi'.")
        sigma_obs_vals.append(sigma_obs)

    sigma_obs_max = sigma_obs_vals[max_index]

    T_fit = np.array(T_values[start_idx:end_idx])
    obs_fit = np.array(obs_values[start_idx:end_idx])
    sigma_obs_vals = np.array(sigma_obs_vals[start_idx:end_idx])

    coeffs, cov = np.polyfit(T_fit, obs_fit, 2, cov=True)
    a, b, c = coeffs
    poly = np.poly1d(coeffs)

    T_c = -coeffs[1] / (2 * coeffs[0]) # -b/2a
    observable_max = poly(T_c) 
    
    var_a = cov[0,0]
    var_b = cov[1,1]
    cov_ab = cov[0,1]

    dTda = b / (2.0 * a**2)
    dTdb = -1.0 / (2.0 * a)
    sigma_Tc = np.sqrt(dTda**2 * var_a + dTdb**2 * var_b + 2 * dTda * dTdb * cov_ab)


    # Plotting
    if plot:
        left = abs(T_c - min(T_fit))
        right = abs(T_c - max(T_fit))
        limits = max(left, right)
        plt.figure()
        T_plot = np.linspace(min(T_fit) - abs(min(T_fit)*0.001), max(T_fit) + abs(max(T_fit)*0.001), 100)
        T_plot = np.linspace(T_c - limits - 0.1*limits, T_c + limits + 0.1*limits, 100)
        plt.plot(T_plot, poly(T_plot), label='Fitted curve')
        plt.errorbar(T_fit, obs_fit, yerr=sigma_obs_vals, fmt="o", color='orange', ecolor="g", capsize=3,  label='Data + errorbars')
        plt.axvline(T_c, color='r', linestyle='--', label=r'$T_c$'+f' = {T_c:.3f}')
        plt.xlabel('Temperature T')
        if observable == "Cv":
            plt.ylabel(r'Heat Capacity $c_v$')
        else:
            plt.ylabel(r'Magnetic Susceptibility $\chi_s$')
        plt.legend()
        plt.tight_layout()
        plt.grid()
        plt.savefig(fig_dir / f'{observable}_Tc_fit_L{L}.pdf')
        plt.close()

    return T_c, observable_max, sigma_obs_max, sigma_Tc

def plot_together(ax_left, ax_right, data, L, observable="Cv"):
    T_values, obs_values = [], []
    for T_str, values in data[str(L)].items():
        T_values.append(float(T_str))
        obs_values.append(values[observable])

    if L > 40:
        linestyle = '--'
    else:
        linestyle = '-'
    # left plot (full range)
    if observable == "chi":
        ax_left.set_yscale('log')
    ax_left.plot(T_values, obs_values, linestyle=linestyle, label=f'{L}', alpha=0.8)
    ax_left.grid(True)

    # right plot (zoomed in at [2.1, 2.4])
    if observable == "chi":
        ax_right.set_yscale('log')
    ax_right.plot(T_values, obs_values, linestyle=linestyle, alpha=0.8)
    ax_right.set_xlim(2.1, 2.4)
    ax_right.grid(True)

def plot_cv_logL(cv_max, L_list):
    """
    Plot Cv_max vs log(L)
    """
    log_L = np.log(L_list)

    plt.figure()
    plt.plot(log_L, cv_max, 'o-')
    plt.xlabel(r'$\ln(L)$')
    plt.ylabel(r'$c_{v_{max}}$')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / 'cv_max_vs_logL.pdf')
    plt.close()
    return None 

# Load data and plot for L
raw_data = load_JSON(json_path)
data = sort_data(raw_data, min_sweeps=min_sweeps)
data_fit = sort_data(raw_data, min_sweeps=min_sweeps_fit)

Tc_Cv_vals = []
Tc_chi_vals = []
sigma_TcCv_vals = []
sigma_TcChi_vals = []
sigma_Cv_vals = []
sigma_chi_vals = []
c_v_max_list = []

for l in L:
    if plot_all:
        plot_Cv_vs_T(data, l)
        plot_chi_vs_T(data, l)
        plot_eps_vs_T(data, l)
        plot_m_vs_T(data, l)
    Tc_Cv, Cv_max, sigma_Cv, sigma_TcCv = Tc_regress(data_fit, l, observable="Cv", plot=plot_all)
    Tc_chi, chi_max, sigma_chi, sigma_TcChi = Tc_regress(data_fit, l, observable="chi", plot=plot_all)
    Tc_Cv_vals.append(Tc_Cv)
    Tc_chi_vals.append(Tc_chi)
    sigma_TcCv_vals.append(sigma_TcCv)
    sigma_TcChi_vals.append(sigma_TcChi)
    sigma_Cv_vals.append(sigma_Cv)
    sigma_chi_vals.append(sigma_chi)
    c_v_max_list.append(Cv_max)
    print(f"L={l}: Tc from Cv = {Tc_Cv:.8f}+-{sigma_TcCv:.4f} (max Cv = {Cv_max:.4f}+-{sigma_Cv:.4f}), Tc from chi = {Tc_chi:.4f}+-{sigma_TcChi:.4f} (max chi = {chi_max:.4f}+-{sigma_chi:.4f})")

plot_cv_logL(c_v_max_list, L)

# 1/L vs Tc plot
Tc_Cv_vals = np.array(Tc_Cv_vals)
Tc_chi_vals = np.array(Tc_chi_vals)
sigma_TcCv_vals = np.array(sigma_TcCv_vals)
sigma_TcChi_vals = np.array(sigma_TcChi_vals)

Tc_vals = 0.5 * (Tc_Cv_vals + Tc_chi_vals)
Tc_errs = 0.5 * np.sqrt(sigma_TcCv_vals**2 + sigma_TcChi_vals**2)

sigma_y = np.sqrt(np.mean(Tc_errs**2))

inv_L = 1.0 / L

# linregress
slope, intercept, r_value, p_val, std_err = linregress(inv_L, Tc_vals)

n = len(L)
inv_L_mean = np.mean(inv_L)
intercept_err = sigma_y * np.sqrt(np.sum(inv_L**2) / (n*np.sum((inv_L - inv_L_mean)**2)))

print(f"Tc(inf) = {intercept:.12f} ± {intercept_err:.5f}")
print(f"The slope m = {slope:.3f} ± {std_err:.3f}")
tc_inf = 2.269185311421
print("According to Lars Onsage's value")
print(f"Tc(inf) ≈ {tc_inf}")

plt.figure()
plt.plot(inv_L, Tc_vals, 'o', label=r'$\hat{T}_c$')
plt.plot(inv_L, intercept + slope * inv_L, 'r--', label=r'$\hat{T}_c$='+f"{intercept:.3f} + {slope:.3f}/L")
plt.xlabel('1/L')
plt.ylabel(r'Critical Temperature $\hat{T}_c$')
plt.legend()
plt.tight_layout()
plt.grid()
plt.savefig(fig_dir / 'Tc_vs_invL.pdf')
plt.close()

data = sort_data(raw_data, min_sweeps=min_sweeps)

# plot-together for multiple observables with stacked legends (Courtesy of ChatGPT to make it look nice)
for observable in ["Cv", "chi", "avg_eps", "avg_mabs"]:
    fig, axes = plt.subplots(ncols=2, figsize=(8, 4.4))
    left, right = axes[0], axes[1]

    for l in L:
        plot_together(left, right, data, l, observable)

    fig.supxlabel('Temperature T', y=0.06)
    if observable == "Cv":
        fig.supylabel(r'$c_V$')
    elif observable == "chi":
        fig.supylabel(r'$\chi_s$')
    elif observable == "avg_eps":
        fig.supylabel(r'$\langle \epsilon \rangle$')
    elif observable == "avg_mabs":
        fig.supylabel(r'$\langle |m| \rangle$')

    # gather unique handles/labels from left axis only (one per L)
    handles, labels = left.get_legend_handles_labels()

    # auto split across top and bottom if many entries
    max_cols = 7  # number of columns per row
    if len(labels) > max_cols:
        mid = int(math.ceil(len(labels) / 2))
        top_h, top_l = handles[:mid], labels[:mid]
        bot_h, bot_l = handles[mid:], labels[mid:]

        leg_top = fig.legend(
            top_h, top_l, loc='upper center', ncol=min(max_cols, len(top_l)),
            bbox_to_anchor=(0.5, 1.04), handlelength=1.0, handletextpad=0.4, columnspacing=0.8,
            frameon=False
        )
        leg_bot = fig.legend(
            bot_h, bot_l, loc='lower center', ncol=min(max_cols, len(bot_l)),
            bbox_to_anchor=(0.5, -0.04), handlelength=1.0, handletextpad=0.4, columnspacing=0.8,
            frameon=False
        )
        fig.add_artist(leg_top)
        fig.subplots_adjust(top=0.86, bottom=0.22)
    else:
        fig.legend(
            handles, labels, loc='lower center', ncol=len(labels),
            bbox_to_anchor=(0.5, -0.04), handlelength=1.0, handletextpad=0.4, columnspacing=0.8,
            frameon=False
        )
        fig.subplots_adjust(bottom=0.22)

    for ax in axes:
        ax.grid(True)

    fig.tight_layout()
    fig.savefig(fig_dir / f'{observable}_vs_T_all_L.pdf')
    plt.close(fig)