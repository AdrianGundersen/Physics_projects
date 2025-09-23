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
    'legend.fontsize': 10,
    'figure.dpi': 300,
})

# Sørg for at output-mappen finnes
os.makedirs("output", exist_ok=True)

path = "output/iterations.txt"
data = np.loadtxt(path)

N, tridiag, random = data.T

# Fit på log-log data
slope_tri, intercept_tri = np.polyfit(np.log10(N), np.log10(tridiag), 1)
slope_rand, intercept_rand = np.polyfit(np.log10(N), np.log10(random), 1)

# Prefaktorer
C_tri = 10**intercept_tri
C_rand = 10**intercept_rand

print(f"Tridiagonal: R(N) ≈ {C_tri:.3f} * N^{slope_tri:.3f}")
print(f"Random:      R(N) ≈ {C_rand:.3f} * N^{slope_rand:.3f}")

# Lag tilpassede kurver (rett linje på log-log)
fit_tri = C_tri * N**slope_tri
fit_rand = C_rand * N**slope_rand

# Plot datapunkter
plt.figure()
plt.loglog(N, tridiag, "o", color="r", markersize=4, label="Tridiagonal data")
plt.loglog(N, random, "s", color="b", markersize=4, label="Random dense data")

# Plot fit-linjer
plt.loglog(N, fit_tri, "-", color="r", label=f"Tridiagonal fit (slope={slope_tri:.2f})")
plt.loglog(N, fit_rand, "-", color="b", label=f"Random fit (slope={slope_rand:.2f})")

plt.xlabel("Matrix size N")
plt.ylabel("Number of Jacobi rotations")
plt.title("Scaling of Jacobi method")
plt.legend()

plt.savefig("output/jacobi_scaling_fit.pdf", bbox_inches="tight")
plt.show()
