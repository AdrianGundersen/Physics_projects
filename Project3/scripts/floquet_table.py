import numpy as np

# ---- constants (SI) ----
e   = 1.602176634e-19
amu = 1.66053906660e-27

def penning_frequencies(V0, d, B0, q=+e, m=40.078*amu):
    wz = np.sqrt(2*q*V0/(m*d**2))
    wc = q*B0/m
    if wc**2 <= 2*wz**2:
        # radial modes complex => unbounded in xy even without drive
        return wz, np.nan, np.nan, wc
    disc = np.sqrt(wc**2 - 2*wz**2)
    w_plus  = 0.5*(wc + disc)
    w_minus = 0.5*(wc - disc)
    return wz, w_plus, w_minus, wc

def rhs_state(t, y, V0, d, B0, f, wV, q, m):
    """y = [x,y,z,vx,vy,vz]."""
    x,yc,z, vx,vy,vz = y
    Ef = V0*(1.0 + f*np.cos(wV*t))/d**2
    Ex, Ey, Ez = Ef*x, Ef*yc, -2.0*Ef*z
    # Lorentz force
    ax = (q/m)*(Ex + vy*B0*1.0 - vz*B0*0.0)  # v×B with B=(0,0,B0)
    ay = (q/m)*(Ey - vx*B0*1.0)
    az = (q/m)*Ez
    return np.array([vx, vy, vz, ax, ay, az], dtype=float)

def rk4_step(fun, t, y, h, *args):
    k1 = fun(t,        y,           *args)
    k2 = fun(t + 0.5*h,y + 0.5*h*k1,*args)
    k3 = fun(t + 0.5*h,y + 0.5*h*k2,*args)
    k4 = fun(t +     h,y +     h*k3,*args)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def growth_rate_per_period(times, signal, periods=100):
    """Least-squares slope of ln|signal| over last `periods` drive periods (per period units)."""
    sig = np.abs(signal)
    sig[sig <= 1e-300] = 1e-300
    L = len(times)
    if L < 10: return 0.0
    # use last window
    t0 = times[-1] - periods
    mask = times >= t0
    x = times[mask]               # measured in drive periods already
    y = np.log(sig[mask])
    # linear fit y = a + b x
    A = np.vstack([np.ones_like(x), x]).T
    a,b = np.linalg.lstsq(A, y, rcond=None)[0]
    return b  # units: per period

def bounded_3d_single_particle(
        V0, d, B0, q=+e, m=40.078*amu,
        f=0.5, wV=2*np.pi*1.0e6,        # rad/s
        T_periods=500, ppw=400,         # simulate 500 drive periods, 400 steps/period
        box_radius=50e-3,               # 50 mm safety box (set to your simulation box)
        growth_tol=1e-4                 # per period threshold
    ):
    """
    Returns dict with bounded flag, growth rates, and mode freqs.
    """
    # frequencies for diagnostics
    wz, w_plus, w_minus, wc = penning_frequencies(V0, d, B0, q=q, m=m)

    # integration setup in drive-period units
    T  = 2*np.pi/wV
    h  = T/ppw
    N  = int(T_periods*ppw)
    t  = 0.0
    # nondeterministic but reasonable small ICs (adjust to your test)
    y = np.array([20e-6, 0.0, 20e-6, 0.0, 25e-6, 0.0])  # [m, m, m, m/s, m/s, m/s]

    times = np.empty(N+1)
    Rhist = np.empty(N+1)
    Zhist = np.empty(N+1)

    times[0] = 0.0
    Rhist[0] = np.hypot(y[0], y[1])
    Zhist[0] = np.abs(y[2])

    escaped = False
    for i in range(1, N+1):
        y = rk4_step(rhs_state, t, y, h, V0, d, B0, f, wV, q, m)
        t += h
        R = np.hypot(y[0], y[1])
        Z = abs(y[2])
        times[i] = i/ppw            # time in *drive periods*
        Rhist[i] = R
        Zhist[i] = Z
        if R > box_radius or Z > box_radius:
            escaped = True
            times = times[:i+1]; Rhist = Rhist[:i+1]; Zhist = Zhist[:i+1]
            break

    # growth rates (per period)
    gR = growth_rate_per_period(times, Rhist, periods=min(100, int(T_periods*0.5)))
    gZ = growth_rate_per_period(times, Zhist, periods=min(100, int(T_periods*0.5)))

    bounded = (not escaped) and (gR <= growth_tol) and (gZ <= growth_tol)

    # quick resonance proximity (optional)
    cands = {}
    if np.isfinite(w_plus):
        cands["2ωz"]     = 2*wz
        cands["ω+±ωz"]   = [abs(w_plus-wz), w_plus+wz]
        cands["2ω+"]     = 2*w_plus
        cands["2ω-"]     = 2*w_minus
        cands["ω++ω-"]   = w_plus + w_minus
        cands["|ω+-ω-|"] = abs(w_plus - w_minus)

    near = []
    for name, val in cands.items():
        vals = val if isinstance(val, (list, tuple)) else [val]
        for vv in vals:
            det = abs(wV - vv)/vv if vv>0 else np.inf
            if det < 0.03:  # within 3%
                near.append((name, vv, det))

    return {
        "bounded": bool(bounded),
        "escaped_geometrically": bool(escaped),
        "growth_rate_R_per_period": float(gR),
        "growth_rate_z_per_period": float(gZ),
        "periods_simulated": float(times[-1]),
        "wz": float(wz), "w_plus": float(w_plus), "w_minus": float(w_minus), "wc": float(wc),
        "near_resonances_3pct": near
    }

# ---------- Example usage ----------
if __name__ == "__main__":
    # Example SI trap numbers (edit): V0=25 mV, d=500 µm, B0=1 T
    V0 = 25e-3
    d  = 500e-6
    B0 = 1.0

    # drive: f and ωV (rad/s). Here 1.6 MHz and f=0.7 as discussed
    wV = 2*np.pi*1.6e6
    out = bounded_3d_single_particle(V0, d, B0, f=0.7, wV=wV,
                                     T_periods=600, ppw=400,
                                     box_radius=1e-2, growth_tol=2e-4)
    print(out)
