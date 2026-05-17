"""
Animate the channel395 LES run as a GIF.
Reads all VTK time snapshots, interpolates instantaneous U_x onto a regular
grid at the z-midplane, and saves images/channel395.gif.
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import griddata

import vtk as _vtk
import vtk.util.numpy_support as vtk_np
import imageio.v2 as imageio

IMAGES   = Path("images")
VTK_DIR  = Path("VTK")
OUT_GIF  = IMAGES / "channel395.gif"

# Flow parameters
Re_tau = 395.0
nu     = 2e-5
u_tau  = Re_tau * nu / 1.0

# ── Collect and sort VTK snapshots ───────────────────────────────────────────
vtk_files = sorted(VTK_DIR.glob("openfoam-channel395_*.vtk"),
                   key=lambda p: int(re.search(r"_(\d+)\.vtk$", p.name).group(1)))
print(f"Found {len(vtk_files)} snapshots: t = "
      + ", ".join(re.search(r"_(\d+)\.vtk$", p.name).group(1) for p in vtk_files))

# ── Build interpolation grid once (shared across all frames) ─────────────────
nx_g, ny_g = 400, 200
xi = np.linspace(0, 4, nx_g)
t_g = np.linspace(0, np.pi, ny_g)
yi  = 1.0 - np.cos(t_g)          # cosine clustering near both walls
XI, YI = np.meshgrid(xi, yi)

# Determine global colour range from first non-zero frame
def read_vtk(path):
    r = _vtk.vtkUnstructuredGridReader()
    r.SetFileName(str(path))
    r.ReadAllScalarsOn(); r.ReadAllVectorsOn()
    r.Update()
    g = r.GetOutput()
    pts = vtk_np.vtk_to_numpy(g.GetPoints().GetData())
    pd  = g.GetPointData()
    def arr(name):
        a = pd.GetArray(name)
        return vtk_np.vtk_to_numpy(a) if a else None
    return pts, arr("U"), arr("UMean")

# z-midplane from first file
pts0, _, _ = read_vtk(vtk_files[0])
z_unique = np.unique(pts0[:, 2])
z_mid    = z_unique[np.argmin(np.abs(z_unique - 1.0))]
mask_z0  = np.abs(pts0[:, 2] - z_mid) < 1e-4
xy_sl    = np.column_stack([pts0[mask_z0, 0], pts0[mask_z0, 1]])

# Sample a mid-run frame for colour range
pts_s, U_s, _ = read_vtk(vtk_files[len(vtk_files)//2])
mask_s = np.abs(pts_s[:, 2] - z_mid) < 1e-4
Ui_s   = griddata(np.column_stack([pts_s[mask_s, 0], pts_s[mask_s, 1]]),
                  U_s[mask_s, 0] / u_tau, (XI, YI), method="linear")
vmin = max(0, np.nanpercentile(Ui_s, 1))
vmax = np.nanpercentile(Ui_s, 99)
levels = np.linspace(vmin, vmax, 50)

# ── Render frames ─────────────────────────────────────────────────────────────
frame_paths = []
for vtk_file in vtk_files:
    t_str = re.search(r"_(\d+)\.vtk$", vtk_file.name).group(1)
    t_val = int(t_str)

    pts, U_inst, UMean = read_vtk(vtk_file)
    mask_z = np.abs(pts[:, 2] - z_mid) < 1e-4
    xy_frame = np.column_stack([pts[mask_z, 0], pts[mask_z, 1]])

    # Use instantaneous U if available, else fall back to UMean
    field = U_inst if U_inst is not None else UMean
    Ux = griddata(xy_frame, field[mask_z, 0] / u_tau, (XI, YI), method="linear")

    fig, ax = plt.subplots(figsize=(10, 3.5))
    cf = ax.contourf(XI, YI, Ux, levels=levels, cmap="RdYlBu_r", extend="both")
    ax.contour(XI, YI, Ux, levels=levels[::5], colors="k",
               linewidths=0.2, alpha=0.35)
    cb = fig.colorbar(cf, ax=ax, pad=0.015, aspect=25)
    cb.set_label("U_x / u_τ", fontsize=9)

    # Walls
    ax.fill_between([0, 4], 0, -0.04, color="dimgray")
    ax.fill_between([0, 4], 2, 2.04,  color="dimgray")
    ax.set_xlim(0, 4); ax.set_ylim(-0.04, 2.04)
    ax.set_xlabel("x (m) — streamwise", fontsize=10)
    ax.set_ylabel("y (m)", fontsize=10)
    ax.set_title(
        f"Channel Reτ = {Re_tau:.0f} — instantaneous U_x⁺  "
        f"(t = {t_val} s  ≈  {t_val * u_tau:.0f} δ/uτ)",
        fontsize=10)

    frame_path = IMAGES / f"_frame_{t_str:>06}.png"
    fig.tight_layout()
    fig.savefig(frame_path, dpi=120)
    plt.close(fig)
    frame_paths.append(frame_path)
    print(f"  frame t={t_val}")

# ── Stitch into GIF ───────────────────────────────────────────────────────────
frames = [imageio.imread(p) for p in frame_paths]
imageio.mimsave(OUT_GIF, frames, duration=0.15, loop=0)
print(f"\nSaved {OUT_GIF}  ({len(frames)} frames)")

# Clean up temp frames
for p in frame_paths:
    p.unlink()
