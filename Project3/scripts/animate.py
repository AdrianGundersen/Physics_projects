import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ---- Last inn data for begge partikler ----
data1 = np.loadtxt("data/pos_vel_1.txt")
data2 = np.loadtxt("data/pos_vel_2.txt")

t  = data1[:,0]   # antar begge har samme tid
x1 = data1[:,1]
y1 = data1[:,2]

x2 = data2[:,1]
y2 = data2[:,2]

# ---- Sett opp plott ----
fig, ax = plt.subplots(figsize=(6,6))
xmin = min(np.min(x1), np.min(x2)) * 1.1
xmax = max(np.max(x1), np.max(x2)) * 1.1
ymin = min(np.min(y1), np.min(y2)) * 1.1
ymax = max(np.max(y1), np.max(y2)) * 1.1

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xlabel(r"x [$\mu$m]")
ax.set_ylabel(r"y [$\mu$m]")
ax.set_title("Particles in xy-plane")

# Punkter og baner for begge partikler
point1, = ax.plot([], [], 'ro', label="Particle 1")
trail1, = ax.plot([], [], 'r-', lw=1, alpha=0.5)

point2, = ax.plot([], [], 'bo', label="Particle 2")
trail2, = ax.plot([], [], 'b-', lw=1, alpha=0.5)

ax.legend()

# ---- Oppdateringsfunksjon ----
def update(fr):
    frame = fr*1000
    # Partikkel 1
    point1.set_data(x1[frame], y1[frame])
    trail1.set_data(x1[:frame], y1[:frame])
    # Partikkel 2
    point2.set_data(x2[frame], y2[frame])
    trail2.set_data(x2[:frame], y2[:frame])
    return point1, trail1, point2, trail2

# ---- Lag animasjon ----
ani = FuncAnimation(
    fig,
    update,
    frames=len(t),
    interval=10,
    blit=True
)

# ---- Vis animasjonen ----
ani.save("data/plot/two_particles.gif", fps=60, dpi=200)
plt.show()

# ---- (valgfritt) lagre som video ----
