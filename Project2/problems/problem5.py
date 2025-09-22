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

# Sørg for at output-mappen finnes
os.makedirs("output", exist_ok=True)

path = "output/iterations.txt"

data = np.loadtxt(path)

N, tridiag, random = data.T

# Plot på log-log skala for å se skalering tydelig
plt.figure()
plt.loglog(N, tridiag, "o-", markersize=4, label="Tridiagonal matrix")
plt.loglog(N, random, "s-", markersize=4, label="Random dense matrix")

plt.xlabel("Matrix size N")
plt.ylabel("Number of Jacobi rotations")
plt.title("Scaling of Jacobi method")
plt.legend()

plt.savefig("output/jacobi_scaling.png", bbox_inches="tight")
plt.show()