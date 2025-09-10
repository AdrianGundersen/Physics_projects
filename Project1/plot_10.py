import numpy as np
import matplotlib.pyplot as plt
import os

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

folder = "output"

data_file = "time_optimized_algo.txt"
data = np.loadtxt(os.path.join(folder, data_file))

opt = data[:,0]
old = data[:,1]

prosent = np.zeros(len(opt))

for i in range(len(opt)):
    prosent[i] = 100 - opt[i] / old[i] * 100
    print(f"{prosent[i]:.4f} %")