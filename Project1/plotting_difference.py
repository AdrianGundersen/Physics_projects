# plotting_difference.py
# Plots the absolute and relative difference between numerical and analytical solutions from a data file.

import os
import numpy as np
import matplotlib.pyplot as plt

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

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir)

data_file = input("Enter the data file name: ")
title = input("Enter the plot title: ")

data = np.loadtxt(os.path.join(file_path, data_file))
