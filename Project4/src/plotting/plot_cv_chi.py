# src/plotting/plot_Cv.py
import numpy as np
import matplotlib as mpl
mpl.use("Agg") # to avoid wayland issues
import matplotlib.pyplot as plt
import json

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

# Load data and plot for L
L = [20]
json_path = ROOT / "test.json"
data = load_JSON(json_path)
for l in L:
    plot_Cv_vs_T(data, l)
    plot_chi_vs_T(data, l)

def find_extremum(data, L, observable="Cv"):
    return None