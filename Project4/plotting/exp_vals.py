# plotting/exp_vals.py
# quick script to calculate average magnetization and energy per spin from data file (used primarily for 2x2 lattice)
import numpy as np
import matplotlib.pyplot as plt
import sys


def exp_vals(eps, mabs, L=2, T=1.0):
    N = L**2
    total_mabs = np.sum(mabs)
    total_eps = np.sum(eps)

    avg_eps = total_eps / len(eps)
    avg_mabs = total_mabs / len(mabs)


    avg_eps2 = np.sum(eps**2) / len(eps)
    avg_mabs2 = np.sum(mabs**2) / len(mabs)

    avg_eps_var = avg_eps2 - avg_eps**2
    avg_mabs_var = avg_mabs2 - avg_mabs**2

    cv = N/(T**2) * avg_eps_var
    chi = N/(T) * avg_mabs_var


    return avg_mabs, avg_eps, cv, chi

def analytical_vals(L=2, T=1.0):
    N = L**2
    k_b = 1.0
    J = 1.0

    beta = 1/(k_b*T)
    Z = 12 + 4*np.cosh(8*J*beta)
    exp_eps_analytical = -32*J/(N*Z) * np.sinh(8*beta*J)
    exp_eps2_analytical = 256*J**2/(N**2*Z) * np.cosh(8*beta*J)
    exp_mabs_analytical = 8/(N*Z) * (2 + np.exp(8*beta*J))
    exp_mabs2_analytical = 32/(N**2*Z) * (1 + np.exp(8*beta*J))
    c_v_analytical = N/(k_b*T**2) * (exp_eps2_analytical - exp_eps_analytical**2)
    chi_analytical = N/(k_b*T) * (exp_mabs2_analytical - exp_mabs_analytical**2)

    return exp_mabs_analytical, exp_eps_analytical, c_v_analytical, chi_analytical


def rel_error(simulated, analytical):
    return np.abs((simulated - analytical) / analytical) * 100  # in percentage

def main(file_path, L, T):
    # Load data from the specified file
    data = np.loadtxt(file_path, delimiter=',', skiprows=0)
    eps, mabs = data[:, 0], data[:, 1]
    
    avg_mabs, avg_eps, cv, chi = exp_vals(eps, mabs, L, T)
    exp_mabs_analytical, exp_eps_analytical, c_v_analytical, chi_analytical = analytical_vals(L, T)

    rel_error_mabs = rel_error(avg_mabs, exp_mabs_analytical)
    rel_error_eps = rel_error(avg_eps, exp_eps_analytical)
    rel_error_cv = rel_error(cv, c_v_analytical)
    rel_error_chi = rel_error(chi, chi_analytical)

    print(f"For file {file_path}:")
    print(f"Average |m|: {avg_mabs}")
    print(f"Average ε: {avg_eps}")
    print(f"Specific Heat Cv: {cv}")
    print(f"Magnetic Susceptibility χ: {chi}")
    
    print("\nAnalytical Values:")
    print(f"Analytical Average |m|: {exp_mabs_analytical}")
    print(f"Analytical Average ε: {exp_eps_analytical}")
    print(f"Analytical Specific Heat Cv: {c_v_analytical}")
    print(f"Analytical Magnetic Susceptibility χ: {chi_analytical}")

    print("\n Relative Errors:")
    print(f"Relative Error |m|: {rel_error_mabs} %")
    print(f"Relative Error ε: {rel_error_eps} %")
    print(f"Relative Error Cv: {rel_error_cv} %")
    print(f"Relative Error χ: {rel_error_chi} %") 

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python exp_vals.py <file_path> <L> <T>")
    else:
        main(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))