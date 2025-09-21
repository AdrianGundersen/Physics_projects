
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

# -----------------------------
# Filstier
# -----------------------------
eigvals_num_10 = "output/problem6_eigenvalues10.csv"
eigvecs_num_10 = "output/problem6_eigenvectors10.csv"
eigvecs_an_10  = "output/problem6_eigenvectors_analytical10.csv"

eigvals_num_100 = "output/problem6_eigenvalues100.csv"
eigvecs_num_100 = "output/problem6_eigenvectors100.csv"
eigvecs_an_100  = "output/problem6_eigenvectors_analytical100.csv"

# -----------------------------
# Last inn data
# -----------------------------
vals_num_10 = np.loadtxt(eigvals_num_10, delimiter=",")
vecs_num_10 = np.loadtxt(eigvecs_num_10, delimiter=",")
vecs_an_10  = np.loadtxt(eigvecs_an_10, delimiter=",")

vals_num_100 = np.loadtxt(eigvals_num_100, delimiter=",")
vecs_num_100 = np.loadtxt(eigvecs_num_100, delimiter=",")
vecs_an_100  = np.loadtxt(eigvecs_an_100, delimiter=",")

def add_boundaries_zero(vecs):
    """Legg til null først og sist langs posisjonsaksen"""
    n, m = vecs.shape  # n = antall punkter, m = antall egenvektorer
    out = np.zeros((n+2, m))
    out[1:-1, :] = vecs
    return out

# Bruk funksjonen på alle egenvektor-matriser
vecs_num_10  = add_boundaries_zero(vecs_num_10)
vecs_an_10   = add_boundaries_zero(vecs_an_10)
vecs_num_100 = add_boundaries_zero(vecs_num_100)
vecs_an_100  = add_boundaries_zero(vecs_an_100)

# X-akse
x_10   = np.linspace(0, 1, vecs_num_10.shape[0])
x_100  = np.linspace(0, 1, vecs_num_100.shape[0])

# -----------------------------
# Plot 1: Alle egenvektorer side-om-side (som før)
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for i in range(vecs_num_10.shape[1]):
    axes[0].plot(x_10, vecs_num_10[:, i], marker="o", label=f"Eigenvector {i+1}")
axes[0].set_title("n=10 eigenvectors")
axes[0].set_xlabel(r"$\hat{x}$")
axes[0].set_ylabel(r"$\vec{v}$")
axes[0].grid(True)
axes[0].legend()

for i in range(vecs_num_100.shape[1]):
    axes[1].plot(x_100, vecs_num_100[:, i], marker="o", label=f"Eigenvector {i+1}")
axes[1].set_title("n=100 eigenvectors")
axes[1].set_xlabel(r"$\hat{x}$")
axes[1].grid(True)
axes[1].legend()

fig.suptitle("Comparison of Eigenvectors")
fig.tight_layout()
fig.savefig("output/eigenvectors_comparison.pdf", bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 2: Egenverdier side-om-side (som før)
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

axes[0].stem(np.arange(1, len(vals_num_10)+1), vals_num_10, basefmt=" ")
axes[0].set_title("n=10 eigenvalues")
axes[0].set_xlabel("Index")
axes[0].set_ylabel("Eigenvalue")
axes[0].grid(True)

axes[1].stem(np.arange(1, len(vals_num_100)+1), vals_num_100, basefmt=" ")
axes[1].set_title("n=100 eigenvalues")
axes[1].set_xlabel("Index")
axes[1].grid(True)

fig.suptitle("Comparison of Eigenvalues")
fig.tight_layout()
fig.savefig("output/eigenvalues_comparison.pdf", bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 3: Numerisk vs analytisk, N=10
# -----------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)

for k in range(vecs_num_10.shape[1]):
    axes[k].plot(x_10, vecs_num_10[:, k], "o-", label="Numerical")
    axes[k].plot(x_10, vecs_an_10[:, k], "--", label="Analytical")
    axes[k].set_title(f"n=10, eigenvector {k+1}")
    axes[k].set_xlabel(r"$\hat{x}$")
    axes[k].grid(True)
    if k == 0:
        axes[k].set_ylabel(r"$v$")
        axes[k].legend()

fig.tight_layout()
fig.savefig("output/eigenvectors_compare_N10.pdf", bbox_inches="tight")
plt.show()

# -----------------------------
# Plot 4: Numerisk vs analytisk, N=100
# -----------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)

for k in range(vecs_num_100.shape[1]):
    axes[k].plot(x_100, vecs_num_100[:, k], label="Numerical")
    axes[k].plot(x_100, vecs_an_100[:, k], "--", label="Analytical")
    axes[k].set_title(f"n=100, eigenvector {k+1}")
    axes[k].set_xlabel(r"$\hat{x}$")
    axes[k].grid(True)
    if k == 0:
        axes[k].set_ylabel(r"$v$")
        axes[k].legend()

fig.tight_layout()
fig.savefig("output/eigenvectors_compare_N100.pdf", bbox_inches="tight")
plt.show()
