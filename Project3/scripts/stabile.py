import numpy as np

q = 1.0
m = 40.078 # Ca+
z0 = 20.0
d = 500.0
V0 = 25.0e-3 * 9.64852558e7
B = 9.64852558e1
omega_z_2 = 2.0 * q * V0 / (m * d**2)
omega_v_2 = 2.11    **2

k = omega_z_2/omega_v_2

def m(k):
    if k < 0:
        return 2 * np.sqrt(k * (k - 1) * (k - 4) / (3 * k - 8))
    
    elif 0 < k < 1/4:
        return 0.25 * (np.sqrt((9 - 4*k)*(13 - 20*k)) - (9 - 4*k))
    
    elif 1/4 <= k < 13/20:
        term = np.sqrt((9 - 4*k)*(13 - 20*k))
        return 0.25 * (9 - 4*k - term), 0.25 * (9 - 4*k + term)  # gir begge ±-løsningene
    
    elif 13/20 <= k < 1:
        return np.sqrt(2 * (k - 1) * (k - 4) * (k - 9) / (k - 5))
    
    elif k > 1:
        return 2 * np.sqrt(k * (k - 1) * (k - 4) / (3*k - 8))
    
    else:
        return np.nan  # for udefinerte punkter

print(f"K    = {k}")
print(f"m(k) = {m(k)}")