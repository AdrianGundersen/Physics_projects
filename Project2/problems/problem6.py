import numpy as np
import matplotlib.pyplot as plt
import os

# Sørg for at output-mappen finnes
os.makedirs("output", exist_ok=True)

# Filstier (to sett med data)
eigvals_file_11 = "output/problem6_eigenvalues.csv"
eigvecs_file_11 = "output/problem6_eigenvectors.csv"

eigvals_file_100 = "output/problem6_eigenvalues100.csv"
eigvecs_file_100 = "output/problem6_eigenvectors100.csv"

# Last inn data
eigenvalues_11 = np.loadtxt(eigvals_file_11, delimiter=",")
eigenvectors_11 = np.loadtxt(eigvecs_file_11, delimiter=",")

eigenvalues_100 = np.loadtxt(eigvals_file_100, delimiter=",")
eigenvectors_100 = np.loadtxt(eigvecs_file_100, delimiter=",")
eigenvectors_100 = eigenvectors_100

# Setting banderies to 0
eigenvectors_11[0, :] = 0     
eigenvectors_11[-1, :] = 0    
eigenvectors_100[0, :] = 0     
eigenvectors_100[-1, :] = 0    

# X-akse
x_11 = np.linspace(0, 1, len(eigenvalues_11))
x_100 = np.linspace(0, 1, len(eigenvalues_100))

# Lag subplots side-om-side for egenvektorer
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for i in range(eigenvectors_11.shape[1]):
    axes[0].plot(x_11, eigenvectors_11[:, i], marker="o", label=f"Eigenvector {i+1}")
axes[0].set_title("n=10 eigenvectors")
axes[0].set_xlabel(r"$\hat{x}$")
axes[0].set_ylabel(r"$\vec{v}$")
axes[0].grid(True)
axes[0].legend()

for i in range(eigenvectors_100.shape[1]):
    axes[1].plot(x_100, eigenvectors_100[:, i], marker="o", label=f"Eigenvector {i+1}")
axes[1].set_title("n=100 eigenvectors")
axes[1].set_xlabel(r"$\hat{x}$")
axes[1].grid(True)
axes[1].legend()

fig.suptitle("Comparison of Eigenvectors")
fig.tight_layout()
fig.savefig("output/eigenvectors_comparison.pdf", bbox_inches="tight")
plt.show()

# Lag subplots side-om-side for egenverdier
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

axes[0].stem(np.arange(1, len(eigenvalues_11)+1), eigenvalues_11, basefmt=" ")
axes[0].set_title("N=11 eigenvalues")
axes[0].set_xlabel("Index")
axes[0].set_ylabel("Eigenvalue")
axes[0].grid(True)

axes[1].stem(np.arange(1, len(eigenvalues_100)+1), eigenvalues_100, basefmt=" ")
axes[1].set_title("N=100 eigenvalues")
axes[1].set_xlabel("Index")
axes[1].grid(True)

fig.suptitle("Comparison of Eigenvalues")
fig.tight_layout()
fig.savefig("output/eigenvalues_comparison.pdf", bbox_inches="tight")
plt.show()
