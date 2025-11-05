import numpy as np
import matplotlib as mpl
mpl.use("Agg") # to avoid wayland issues
import matplotlib.pyplot as plt
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

fig_dir = ROOT / "data/figures"
fig_dir.mkdir(parents=True, exist_ok=True)

filepath = ROOT / "data/output/critical_temperature.txt"

unsorted_data = np.loadtxt(filepath, delimiter=',')
data = unsorted_data[unsorted_data[:,0].argsort()]
L = data[:,0].astype(int)
T_C_V = data[:,1]
T_SUSC = data[:,3]

T_L = np.zeros(len(L))
for i in range(len(L)):
    T_L[i] = (T_C_V[i] + T_SUSC[i]) / 2

x = 1/np.array(L)
y = T_L
print(y)

a, T_c_infty  = np.polyfit(x, y, 1)

print(f"Extrapolated critical temperature T_c(inf) = {T_c_infty:0.6f}")

x_fit = np.linspace(np.min(x), np.max(x), 100)
y_fit = a * x_fit + T_c_infty

plt.figure()          
plt.plot(x, y, 'o-', label=r'$T_c$ points')
plt.plot(x_fit, y_fit, 'r--', label=r'$T_c(L)$')
plt.xlabel(r'$1/L$')
plt.ylabel(r'$T_c(L)$')
plt.title(r'$T_c(L=\infty) = %.6f$' % T_c_infty)
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "T_c_finite_size_scaling.pdf")
plt.close()

plt.figure()          
plt.plot(L, y, 'o-', label=r'$T_c$ points')
plt.xlabel(r'$L$')
plt.ylabel(r'$T_c(L)$')
plt.title('Critical temperature vs Lattice Size')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "critical_temperature_vs_lattice_size.pdf")
plt.close()