import numpy as np
import matplotlib as mpl
mpl.use("Agg") # to avoid wayland issues
import matplotlib.pyplot as plt
from pathlib import Path
import math
import json


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

N = 4
k_b = 1.0
T = "2.000000" #change depending on the file you wish to analyze
J = 1.0
p10 = 1000 # amount of points 


# Finds project root
ROOT = Path(__file__).resolve().parents[2]
data_root = ROOT/f"Project4/data/outputs/L2_T={T}_spin=random.txt"
T = float(T)
data = np.loadtxt(data_root, delimiter=',', skiprows=0)
eps, mabs =  data[:,0], data[:,1]
eps2, mabs2 = eps**2, mabs**2


# Load the json file
with open(ROOT / "Project4/data/outputs/tsweep_results.json", "r") as f:   
    data = json.load(f)

fig_dir = ROOT / "Project4/data/figures/L2x2"
fig_dir.mkdir(parents=True, exist_ok=True)




#--------------- Data extraction --------------
L = 2
L2_data = data[str(L)]   

T_values = []
Cv_values = []
chi_values = []
eps_values = []
mabs_values = []
eps2_values = []
mabs2_values = []

for T_str, values in data[str(L)].items():
    T_values.append(float(T_str))
    Cv_values.append(values["Cv"]/N)
    chi_values.append(values["chi"]/N)
    eps_values.append(values["avg_eps"])
    mabs_values.append(values["avg_mabs"])
    eps2_values.append(values["avg_eps2"])
    mabs2_values.append(values["avg_mabs2"])



# p10 = int(math.log10(len(eps)))
exp_eps, exp_eps2, exp_mabs, exp_mabs2 = np.zeros(p10), np.zeros(p10), np.zeros(p10), np.zeros(p10) 
part = int(len(eps)//p10)
for i in range(0,p10):
    exp_eps[i] = np.average(eps[:part*(i+1)])
    exp_eps2[i] = np.average(eps2[:part*(i+1)])
    exp_mabs[i] = np.average(mabs[:part*(i+1)])
    exp_mabs2[i] = np.average(mabs2[:part*(i+1)])

c_v = N/(k_b*T**2) * (exp_eps2 - exp_eps**2)
chi = N/(k_b*T) * (exp_mabs2 - exp_mabs**2)

beta = 1/(k_b*T)
Z = 12 + 4*np.cosh(8*J*beta)
exp_eps_analytical = -32*J/(N*Z) * np.sinh(8*beta*J)
exp_eps2_analytical = 256*J**2/(N**2*Z) * np.cosh(8*beta*J)
exp_mabs_analytical = 8/(N*Z) * (2 + np.exp(8*beta*J))
exp_mabs2_analytical = 32/(N**2*Z) * (1 + np.exp(8*beta*J))
c_v_analytical = N/(k_b*T**2) * (exp_eps2_analytical - exp_eps_analytical**2)
chi_analytical = N/(k_b*T) * (exp_mabs2_analytical - exp_mabs_analytical**2)



n_cycles = np.linspace(1, len(eps), p10)


rel_error_eps = np.abs((exp_eps - exp_eps_analytical)/exp_eps_analytical)
rel_error_mabs = np.abs((exp_mabs - exp_mabs_analytical)/exp_mabs_analytical)
rel_error_cv = np.abs((c_v - c_v_analytical)/c_v_analytical)
rel_error_chi = np.abs((chi - chi_analytical)/chi_analytical)

#-----------------------------------------------------------



#----------------- Plotting ----------------------------
# Energy per spin relative error evolution
plt.figure()
plt.plot(n_cycles, rel_error_eps, label='relative error')
plt.xlabel(f'Cycles for T = {T}')
plt.ylabel(r'Relative error $\langle \epsilon \rangle$')
plt.grid(True, alpha=0.6)
plt.legend()
plt.savefig(fig_dir/f"eps_evolution_T={T}.pdf", bbox_inches='tight')
plt.close()


# Magnetization per spin relative error evolution
plt.figure()
plt.plot(n_cycles, rel_error_mabs, label='relative error')
plt.xlabel(f'Cycles for T = {T}')
plt.ylabel(r'Relative error $\langle |m| \rangle$')
plt.grid(True, alpha=0.6)
plt.legend()
plt.savefig(fig_dir/f"mabs_evolution_T={T}.pdf", bbox_inches='tight')
plt.close()      

# Heat capacity per spin relative error evolution
plt.figure()
plt.plot(n_cycles, rel_error_cv, label='relative error')
plt.xlabel(f'Cycles for T = {T}')
plt.ylabel(r'Relative error $C_V/N$')
plt.grid(True, alpha=0.6)
plt.legend()
plt.savefig(fig_dir/f"cv_evolution_T={T}.pdf", bbox_inches='tight')
plt.close()         


# Susceptibility per spin relative error evolution
plt.figure()
plt.plot(n_cycles, rel_error_chi, label ='relative error')
plt.grid(True, alpha=0.6)
plt.xlabel(f'Cycles for T = {T}')
plt.ylabel(r'Relative error $\chi/N$')
plt.legend()
plt.savefig(fig_dir/f"chi_evolution_T={T}.pdf", bbox_inches='tight')
plt.close()


#all in same plot
plt.figure()
plt.plot(n_cycles, (rel_error_eps), label=r'$\langle \epsilon \rangle$')
plt.plot(n_cycles, (rel_error_mabs), label=r'$\langle |m| \rangle$')
plt.plot(n_cycles, (rel_error_cv), label=r'$C_V/N$')
plt.plot(n_cycles, (rel_error_chi), label=r'$\chi/N$')
plt.xscale('log')   
# plt.yscale('log')
plt.xlabel(f'Cycles for T = {T}')
plt.ylabel('Relative error')
plt.grid(True, alpha=0.6)
plt.legend()
plt.savefig(fig_dir/f"all_rel_error_evolution_T={T}.pdf", bbox_inches='tight')
plt.close()




#-------------- T range -----------------
#--------------- Data extraction --------------
T_values = np.array(T_values)
Cv_values = np.array(Cv_values)
chi_values = np.array(chi_values)
eps_values = np.array(eps_values)
mabs_values = np.array(mabs_values) 


beta = 1/(k_b*T_values)
Z = 12 + 4*np.cosh(8*J*beta)
exp_eps_analytical = -32*J/(N*Z) * np.sinh(8*beta*J)
exp_eps2_analytical = 256*J**2/(N**2*Z) * np.cosh(8*beta*J)
exp_mabs_analytical = 8/(N*Z) * (2 + np.exp(8*beta*J))
exp_mabs2_analytical = 32/(N**2*Z) * (1 + np.exp(8*beta*J))
c_v_analytical = N/(k_b*T_values**2) * (exp_eps2_analytical - exp_eps_analytical**2)
chi_analytical = N/(k_b*T_values) * (exp_mabs2_analytical - exp_mabs_analytical**2)


rel_error_eps = np.abs((eps_values - exp_eps_analytical)/exp_eps_analytical)
rel_error_mabs = np.abs((mabs_values - exp_mabs_analytical)/exp_mabs_analytical)
rel_error_cv = np.abs((Cv_values - c_v_analytical/N)/ (c_v_analytical/N))
rel_error_chi = np.abs((chi_values - chi_analytical/N)/(chi_analytical/N))
    
   
#------------------------------------------------------------------

#eps vs T relative error
plt.figure()
plt.plot(T_values, rel_error_eps, marker='o', linestyle='-', label='relative error')
plt.grid(True, alpha=0.6)
plt.xlabel(f'Temperature T')
plt.ylabel(r'$\langle \epsilon \rangle$')
plt.legend()
plt.savefig(fig_dir/f"eps_vs_T_rel_error_L=2.pdf", bbox_inches='tight')
plt.close()

#mabs vs T relative error
plt.figure()
plt.plot(T_values, rel_error_mabs, marker='o', linestyle='-', label='relative error')
plt.grid(True, alpha=0.6)
plt.xlabel(f'Temperature T')
plt.ylabel(r'$\langle |m| \rangle$')
plt.grid(True, alpha=0.6)
plt.legend()
plt.savefig(fig_dir/f"mabs_vs_T_rel_error_L=2.pdf", bbox_inches='tight')
plt.close() 

#Cv vs T relative error
plt.figure()
plt.plot(T_values, rel_error_cv, marker='o', linestyle='-', label='relative error')
plt.grid(True, alpha=0.6)
plt.xlabel(f'Temperature T')
plt.ylabel(r'$C_V/N$')
plt.legend()
plt.savefig(fig_dir/f"cv_vs_T_rel_error_L=2.pdf", bbox_inches='tight')
plt.close() 

#chi vs T relative error
plt.figure()
plt.plot(T_values, rel_error_chi, marker='o', linestyle='-', label ='relative error')
plt.grid(True, alpha=0.6)
plt.xlabel(f'Temperature T')
plt.ylabel(r'$\chi/N$')
plt.legend()
plt.savefig(fig_dir/f"chi_vs_T_rel_error_L=2.pdf", bbox_inches='tight')
plt.close()

#all relative errors in same plot
plt.figure()
plt.plot(T_values, (rel_error_eps), marker='o', linestyle='-', label=r'$\langle \epsilon \rangle$')
plt.plot(T_values, (rel_error_mabs), marker='o', linestyle='-', label=r'$\langle |m| \rangle$')
plt.plot(T_values, (rel_error_cv), marker='o', linestyle='-', label=r'$C_V/N$')
plt.plot(T_values, (rel_error_chi), marker='o', linestyle='-', label=r'$\chi/N$')
plt.xscale('log')
# plt.yscale('log')
plt.xlabel(f'Temperature T')
plt.ylabel('Relative error')
plt.grid(True, alpha=0.6)
plt.legend()
plt.savefig(fig_dir/f"all_rel_error_vs_T_L=2.pdf", bbox_inches='tight')
plt.close() 

# Cv vs T
plt.figure()
plt.plot(T_values, Cv_values, linestyle='-', label='Numbercal')
plt.plot(T_values, c_v_analytical/N, linestyle='--', label='Analytical')
plt.xlabel('Temperature T')
plt.ylabel('Heat Capacity Cv')
plt.legend()
plt.tight_layout()
plt.grid(True, alpha=0.6)
plt.savefig(fig_dir / f'Cv_vs_T_L=2.pdf')
plt.close()

# Chi vs T
plt.figure()
plt.plot(T_values, chi_values, linestyle='-', label='Numerical')
plt.plot(T_values, chi_analytical/N, linestyle='--', label='Analytical')
plt.xlabel(r'Temperature $T$')
plt.ylabel(r'Susceptibility $\chi$')
plt.legend()
plt.tight_layout()
plt.grid(True, alpha=0.6)
plt.savefig(fig_dir / f'chi_over_N_vs_T_L=2.pdf')
plt.close()

# eps vs T
plt.figure()
plt.plot(T_values, eps_values, linestyle='-', label='Numerical')
plt.plot(T_values, exp_eps_analytical, linestyle='--', label='Analytical')
plt.xlabel(r'Temperature $T$')
plt.ylabel(r'Energy per spin $\langle \epsilon \rangle$')
plt.legend()
plt.tight_layout()
plt.grid(True, alpha=0.6)
plt.savefig(fig_dir / f'eps_vs_T_L=2.pdf')
plt.close()

# mabs vs T
plt.figure()
plt.plot(T_values, mabs_values, linestyle='-', label='Numerical')
plt.plot(T_values, exp_mabs_analytical, linestyle='--', label='Analytical')
plt.xlabel(r'Temperature $T$')
plt.ylabel(r'Absolute Magnetization per spin $\langle |m| \rangle$')
plt.legend()
plt.tight_layout()
plt.grid(True, alpha=0.6)
plt.savefig(fig_dir / f'mabs_vs_T_L=2.pdf')
plt.close()

