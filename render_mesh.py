"""
Render mesh schematic for the channel395 LES case.
Output: images/mesh.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

IMAGES = Path("images")

# Mesh parameters: 40 x 50 x 30 (two blocks of 40x25x30), domain 4 x 2 x 2
# Wall-normal grading: expansion ratio 10.7 (bottom block), 0.093 (top block)
# y-cell heights computed from grading

def graded_y(n, L, r):
    """Return cell-face y-positions for n cells, length L, grading ratio r."""
    if abs(r - 1) < 1e-6:
        return np.linspace(0, L, n + 1)
    # geometric series: first cell * (1 + r + r^2 + ... + r^(n-1)) = L
    dy0 = L * (r - 1) / (r**n - 1)
    faces = [0.0]
    dy = dy0
    for _ in range(n):
        faces.append(faces[-1] + dy)
        dy *= r
    return np.array(faces)

# Bottom block: y=0 to 1, 25 cells, grading 10.7 (cells grow away from wall)
y_bot = graded_y(25, 1.0, 10.7**(1/24))
# Top block: y=1 to 2, 25 cells, grading 0.0934 (cells shrink toward top wall)
y_top = graded_y(25, 1.0, 0.0934**(1/24)) + 1.0

y_all = np.concatenate([y_bot, y_top[1:]])

# Uniform in x (4m, 40 cells) and z (2m, 30 cells)
x_all = np.linspace(0, 4, 41)
z_all = np.linspace(0, 2, 31)

# Wall units: y+ = y * u_tau / nu
u_tau = 0.0079   # m/s
nu    = 2e-5
yp_all = y_all * u_tau / nu

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# ── Left: x-y plane cross-section (every 4th x-line, every y-line) ───────────
ax = axes[0]
for xv in x_all[::4]:
    ax.axvline(xv, color="steelblue", lw=0.5, alpha=0.7)
for yv in y_all:
    ax.axhline(yv, color="steelblue", lw=0.5, alpha=0.7)
ax.axhline(0, color="k", lw=2)
ax.axhline(2, color="k", lw=2)
ax.set_xlim(0, 4)
ax.set_ylim(0, 2)
ax.set_xlabel("x (m)  — streamwise", fontsize=10)
ax.set_ylabel("y (m)  — wall-normal", fontsize=10)
ax.set_title("x-y cross-section  (flow left → right)\n"
             "40×50 cells shown  |  wall-normal grading visible", fontsize=9)
yp1 = yp_all[1]   # first cell height in wall units
ax.annotate(f"bottom wall  y⁺₁ = {yp1:.1f}", xy=(2, y_all[1]), xytext=(2, 0.25),
            ha="center", fontsize=8, color="darkred",
            arrowprops=dict(arrowstyle="->", color="darkred"))
ax.annotate(f"top wall  y⁺₁ = {yp1:.1f}", xy=(2, 2 - y_all[1]), xytext=(2, 1.75),
            ha="center", fontsize=8, color="darkred",
            arrowprops=dict(arrowstyle="->", color="darkred"))
ax.text(0.2, 1.0, "Periodic\nin x and z", ha="left", va="center",
        fontsize=8, color="navy",
        bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", ec="navy", lw=0.5))

# ── Right: near-wall zoom (y+ < 50) ──────────────────────────────────────────
ax2 = axes[1]
yp_zoom = yp_all[yp_all <= 60]
y_zoom  = y_all[:len(yp_zoom)]

for xv in x_all[::4]:
    ax2.axvline(xv, color="steelblue", lw=0.7, alpha=0.8)
for yv in y_zoom:
    ax2.axhline(yv, color="steelblue", lw=0.7, alpha=0.8)
ax2.axhline(0, color="k", lw=2)

# y+ reference lines
for yp_ref, label in [(5, "y⁺=5\n(sublayer)"), (30, "y⁺=30\n(log layer)")]:
    y_ref = yp_ref * nu / u_tau
    ax2.axhline(y_ref, color="darkorange", lw=0.8, ls="--")
    ax2.text(3.9, y_ref + 0.001, label, ha="right", fontsize=7, color="darkorange")

ax2.set_xlim(0, 4)
ax2.set_ylim(0, yp_zoom[-1] * nu / u_tau)
ax2.set_xlabel("x (m)", fontsize=10)
ax2.set_ylabel("y (m)", fontsize=10)
ax2_r = ax2.twinx()
ax2_r.set_ylim(0, yp_zoom[-1])
ax2_r.set_ylabel("y⁺", fontsize=10, color="darkred")
ax2_r.tick_params(colors="darkred")
ax2.set_title(f"Near-wall zoom  (y⁺ < {yp_zoom[-1]:.0f})\n"
              f"First cell y⁺₁ = {yp_all[1]:.2f}  —  well inside viscous sublayer", fontsize=9)

fig.suptitle("Channel Reτ = 395 — LES mesh  (40 × 50 × 30 = 60 000 cells)", fontsize=11)
fig.tight_layout()
fig.savefig(IMAGES / "mesh.png", dpi=150)
plt.close(fig)
print("mesh.png done")
