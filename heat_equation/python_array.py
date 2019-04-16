import numpy as np


def heat(t, D, v, Nx, Nt):
    dt = t / Nt
    dx = 1 / Nx
    u = np.ones(Nx)
    u_next = np.zeros(Nx)
    for it in range(0, Nt - 1):
        u_next[1:Nx-2] = (
              u[1:Nx-2]
              + D * dt / (dx * dx) * (u[2:Nx-1] - 2 * u[1:Nx-2]+u[0:Nx-3])
              - v * dt / dx * (u[2:Nx-1] - u[1:Nx-2])
        )
        u[:] = u_next[:]
    return u


Nx = 1000
Nt = 1000000
D = 1.14
v = .85
t = 0.1
print(heat(t, D, v, Nx, Nt)[2])
