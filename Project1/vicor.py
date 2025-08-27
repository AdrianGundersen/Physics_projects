import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as cons
from scipy.stats import chi2
import pandas as pd

# Oppdaterer mpl parametere
plt.rcParams.update({
    'font.size': 14,
    'figure.figsize': (6, 4),
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2,
    'legend.fontsize': 10,
    'figure.dpi': 300,
})
# Opprett mappen "figure" dersom den ikke eksisterer
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("figure", exist_ok=True)