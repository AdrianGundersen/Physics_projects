import numpy as np

def floquet_table(omega_V, f):
    q = 1.0
    m = 40.078 # Ca+
    z0 = 20.0
    d = 500.0
    V0 = 25.0e-3 * 9.64852558e7
    B = 9.64852558e1
    omega_z = np.sqrt(2.0 * q * V0 / (m * d**2))
    k = omega_z**2 / omega_V**2

    def s1(x): return np.sqrt((9-4*x)*(13-20*x))
    def s2(x): return np.sqrt(2*(x-1)*(x-4)*(x-9)/(x-5))
    def s3(x): return np.sqrt(x*(x-1)*(x-4)/(3*x-8))

    m_minus = None

    if k < 0:
        r = s3(k);  m_minus, m_plus = -2*r,  2*r
    elif 0 < k < 1/4:
        r = s1(k);  m_minus, m_plus = (r - (9 - 4*k))/4, (r + (9 - 4*k))/4
    elif 1/4 <= k < 13/20:
        r = s1(k);  m_minus, m_plus = ((9 - 4*k) - r)/4, ((9 - 4*k) + r)/4
    elif 13/20 <= k < 1:
        r = s2(k);  m_minus, m_plus = -r, r
    else:  # k >= 1
        r = s3(k);  m_minus, m_plus = -2*r, 2*r

    trapped = (f >= m_minus) and (f <= m_plus)

    print(f"{'Trapped' if trapped else 'Not trapped'} for omega_V={omega_V:.6f}, f={f:.3f} "
          f"(k={k:.6f}, m-={m_minus:.6f}, m+={m_plus:.6f})")

    return k, m_minus, m_plus


omega_V_list = np.linspace(0.2,2.5,50) # MHz
f_tuple = (0.3, 0.5, 0.7)


for i in f_tuple:
    print(f"\n--- f={i} ---")
    for omega_V in omega_V_list:
        floquet_table(omega_V, i)
