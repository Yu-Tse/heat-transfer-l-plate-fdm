#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation, PillowWriter

# ─── physical / numerical parameters ─────────────────────────────
dt, Tend = 0.1, 5_000.0                       # s
time     = np.arange(0.0, Tend + dt, dt)
L        = 0.2                                # m
T0, Ts, T_inf = 300.0, 400.0, 300.0           # K
h, q_dot = 20.0, 20_000.0                     # W/m²K , W/m³
dx = dy = 5e-3                                # m
k, rho, Cp = 15.1, 8_055.0, 480.0             # material props

Bi   = h * dx / k
alpha = k / (rho * Cp)
F0    = alpha * dt / dx**2
n     = int(L / dx) + 1                       # 41 nodes / side
size  = n * n                                 # 1681 unknowns

# q̇ source term from the MATLAB derivation
N1 = 2 * dt * q_dot / (rho * Cp * dx)

# ─── build A (dense 1681×1681) and B (1681) ──────────────────────
A = np.zeros((size, size))
B = np.zeros(size)

for j in range(1, size + 1):        # 1-based index to match MATLAB
    col = (j - 1) // 41 + 1         # MATLAB’s “row” = column in x-direction
    mod41 = j % 41 if j % 41 else 41

    # ---------- 17 node types from textbook ----------
    if col == 1:                                   # first column (x = 0)
        if mod41 == 1:                             # node-1 (Dirichlet)
            B[j-1] = Ts
        elif 2 <= mod41 <= 40:                     # node-2
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = F0
            A[j-1, j    ] = F0
            A[j-1, j+40]  = 2*F0
            B[j-1]        = N1
        else:                                      # node-3
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = 2*F0
            A[j-1, j+40]  = 2*F0
            B[j-1]        = N1

    elif 2 <= col <= 20:                           # columns 2…20
        if mod41 == 1:                             # node-4
            B[j-1] = Ts
        elif 2 <= mod41 <= 40:                     # node-5
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = F0
            A[j-1, j    ] = F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0
        else:                                      # node-6
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = 2*F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0

    elif col == 21:                                # column 21
        if mod41 == 1:                             # node-7
            B[j-1] = Ts
        elif 2 <= mod41 <= 20:                     # node-8
            A[j-1, j-1]   = 1 - 4*F0 - 2*Bi*F0
            A[j-1, j-2]   = F0
            A[j-1, j    ] = F0
            A[j-1, j-42]  = 2*F0
            B[j-1]        = 2*Bi*F0*T_inf
        elif mod41 == 21:                          # node-9
            A[j-1, j-1]   = 1 - 4*F0 - (2/3)*Bi*F0
            A[j-1, j-2]   = (2/3)*F0
            A[j-1, j    ] = (4/3)*F0
            A[j-1, j-42]  = (4/3)*F0
            A[j-1, j+40]  = (2/3)*F0
            B[j-1]        = (2/3)*Bi*F0*T_inf
        elif 22 <= mod41 <= 40:                    # node-10
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = F0
            A[j-1, j    ] = F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0
        else:                                      # node-11
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = 2*F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0

    elif 22 <= col <= 40:                          # columns 22…40
        if mod41 == 21:                            # node-12
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j    ] = 2*F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0
        elif 22 <= mod41 <= 40:                    # node-13
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = F0
            A[j-1, j    ] = F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0
        else:                                      # node-14
            A[j-1, j-1]   = 1 - 4*F0
            A[j-1, j-2]   = 2*F0
            A[j-1, j-42]  = F0
            A[j-1, j+40]  = F0

    else:                                          # col == 41
        if mod41 == 21:                            # node-15
            A[j-1, j-1]   = 1 - 4*F0 - 2*Bi*F0
            A[j-1, j    ] = 2*F0
            A[j-1, j-42]  = 2*F0
            B[j-1]        = 2*F0*Bi*T_inf
        elif 22 <= mod41 <= 40:                    # node-16
            A[j-1, j-1]   = 1 - 4*F0 - 2*Bi*F0
            A[j-1, j-2]   = F0
            A[j-1, j    ] = F0
            A[j-1, j-42]  = 2*F0
            B[j-1]        = 2*Bi*F0*T_inf
        else:                                      # node-17
            A[j-1, j-1]   = 1 - 4*F0 - 2*Bi*F0
            A[j-1, j-2]   = 2*F0
            A[j-1, j-42]  = 2*F0
            B[j-1]        = 2*F0*Bi*T_inf

# ─── initial temperature field (41 × 41) ─────────────────────────
T_grid = np.full((41, 41), T0)
T_grid[0:21, 0]      = Ts     # west boundary
T_grid[0:21, 20:41]  = 0.0    # initial value in the cut-out
T_vec = T_grid.ravel(order='F')   # Fortran order, like MATLAB

# ─── time marching until steady state ────────────────────────────
frames, times = [], []
tol           = 1e-4
capture_every = 500            # save every N steps

for step, t in enumerate(time):
    T_new = A @ T_vec + B
    if step % capture_every == 0:
        frames.append(T_new.copy())
        times .append(t)
    if np.all(np.abs(T_new - T_vec) < tol):
        frames.append(T_new.copy())
        times .append(t)
        break
    T_vec = T_new

# ─── helper: embed 1681-vector back to 41×41 and mask cut-out ────
def embed(vec):
    M = vec.reshape(41, 41, order='F')
    M[:21, 20:41] = np.nan     # mask the L-shape cut-out
    return M

full_frames = [embed(v) for v in frames]
steady_Z    = full_frames[-1]

# ─── static 3-D plot (PNG) ───────────────────────────────────────
fig_static = plt.figure(figsize=(6.5, 4.8))
ax_s  = fig_static.add_subplot(111, projection='3d')
X, Y  = np.meshgrid(np.linspace(0, L, 41),
                    np.linspace(0, L, 41))
surf  = ax_s.plot_surface(X, Y, steady_Z, cmap='hot', linewidth=0,
                          antialiased=False)
ax_s.set(xlabel='x [m]', ylabel='y [m]', zlabel='T [K]',
         title='Steady-state (Explicit FDM)')
cbar = fig_static.colorbar(surf, shrink=0.8, aspect=8, pad=0.1)
cbar.set_label('Temperature [K]')
fig_static.tight_layout()
fig_static.savefig('steady_state.png', dpi=300)

# ─── 3-D rotating animation (one full turn) ─────────────────────
fig_rot = plt.figure(figsize=(6.5, 4.8))
ax_r    = fig_rot.add_subplot(111, projection='3d')
surf_r  = ax_r.plot_surface(X, Y, steady_Z, cmap='hot', linewidth=0,
                            antialiased=False)
ax_r.set(xlabel='x [m]', ylabel='y [m]', zlabel='T [K]',
         title='Steady-state – rotating view')
mappable = cm.ScalarMappable(cmap='hot')
mappable.set_array(steady_Z[np.isfinite(steady_Z)])
fig_rot.colorbar(mappable, ax=ax_r, shrink=0.8, aspect=8, pad=0.1)

def rotate(angle):
    """Rotate the view around z-axis."""
    ax_r.view_init(elev=30, azim=angle)
    return surf_r,

ani_rot = FuncAnimation(fig_rot, rotate,
                        frames=np.linspace(0, 360, 181),  # 2° per frame
                        blit=True)
ani_rot.save('steady_rotation.gif',
             writer=PillowWriter(fps=20))

# ─── 2-D temperature evolution GIF ───────────────────────────────
fig2d, ax2d = plt.subplots()
im    = ax2d.imshow(full_frames[0], origin='lower',
                    extent=[0, L, 0, L], cmap='hot')
cbar2 = fig2d.colorbar(im, ax=ax2d, shrink=0.8, aspect=8)
cbar2.set_label('Temperature [K]')
txt   = ax2d.text(0.02, 0.95, '', transform=ax2d.transAxes,
                  color='w', fontsize=10)

def update2d(i):
    im.set_data(full_frames[i])
    txt.set_text(f't = {times[i]:.0f} s')
    return im, txt

ani_2d = FuncAnimation(fig2d, update2d,
                       frames=len(full_frames), interval=100, blit=True)
ani_2d.save('temp_evolution.gif', writer=PillowWriter(fps=10))

print('✔  Outputs generated:',
      'steady_state.png, steady_rotation.gif, temp_evolution.gif')
