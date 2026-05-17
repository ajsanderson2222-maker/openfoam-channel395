"""
Post-processing for the channel395 LES case.
Validates mean velocity and Reynolds stresses against Moser, Kim & Mansour (1999) DNS.
Produces:
  images/convergence.png  — mean streamwise velocity history
  images/validation.png   — U+ vs y+, k+ vs y+, -<u'v'>+ vs y+
  images/contours.png     — x-y slice of mean and instantaneous streamwise velocity
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

IMAGES = Path("images")
IMAGES.mkdir(exist_ok=True)
DNS    = Path("validation/dns")
PP     = Path("postProcessing/layerAverage")

# Flow parameters
Re_tau = 395.0
nu     = 2e-5       # kinematic viscosity [m^2/s]
delta  = 1.0        # channel half-height [m]
u_tau  = Re_tau * nu / delta   # friction velocity = 0.0079 m/s

# ── DNS reference data ────────────────────────────────────────────────────────
# Moser, Kim & Mansour (1999), Re_tau=392.24 (normalised by u_tau and delta)
# Columns: y, y+, Umean, dUmean/dy, Wmean, dWmean/dy, Pmean
dns_means  = np.loadtxt(DNS / "chan395_means.dat",  comments="#")
# Columns: y, y+, R_uu, R_vv, R_ww, R_uv, R_uw, R_vw
dns_stress = np.loadtxt(DNS / "chan395_stress.dat", comments="#")

# DNS is normalised by u_tau and delta; only take first half (y=0 to 1)
dns_yp  = dns_means[:, 1]          # y+
dns_Up  = dns_means[:, 2]          # U+ = U/u_tau
dns_yp_s   = dns_stress[:, 1]
dns_uu  = dns_stress[:, 2]         # <u'u'>+
dns_vv  = dns_stress[:, 3]         # <v'v'>+
dns_ww  = dns_stress[:, 4]         # <w'w'>+
dns_uv  = dns_stress[:, 5]         # <u'v'>+ (negative in channel)
dns_k   = 0.5 * (dns_uu + dns_vv + dns_ww)

# Only plot half-channel (y+ = 0 to Re_tau)
mask = dns_yp <= Re_tau + 5
dns_yp  = dns_yp[mask]
dns_Up  = dns_Up[mask]
mask_s  = dns_yp_s <= Re_tau + 5
dns_yp_s = dns_yp_s[mask_s]
dns_k   = dns_k[mask_s]
dns_uv  = dns_uv[mask_s]

# ── Load OpenFOAM layer-averaged profiles ────────────────────────────────────
# postProcessing/layerAverage/  — written every timestep
# Columns from graphLayerAverage: y, pMean, pPrime2Mean, UMean_x, UMean_y, UMean_z,
#                                  UPrime2Mean_xx, xy, xz, yy, yz, zz, k
pp_dir = PP
times  = sorted(pp_dir.iterdir(), key=lambda p: float(p.name))

# Use the last written time (fully averaged)
latest = times[-1]
# columns: y, pMean, pPrime2Mean, UMean_x, UMean_y, UMean_z,
#          UPrime2Mean_xx, xy, xz, yy, yz, zz, k
d = np.loadtxt(latest / "layerAverage.xy", comments="#")
print(f"Loaded {d.shape[0]} wall-normal points from t={latest.name}")

y_cfd  = d[:, 0]                   # physical y from bottom wall [m]
Ux_cfd = d[:, 3]                   # UMean_x  [m/s]
uu_cfd = d[:, 6]                   # UPrime2Mean_xx
vv_cfd = d[:, 9]                   # UPrime2Mean_yy
ww_cfd = d[:, 11]                  # UPrime2Mean_zz
uv_cfd = d[:, 7]                   # UPrime2Mean_xy
k_cfd  = d[:, 12]                  # k = 0.5*tr(UPrime2Mean)

# Convert to wall units (y+ from bottom wall)
yp_cfd  = y_cfd * u_tau / nu
Up_cfd  = Ux_cfd / u_tau
kp_cfd  = k_cfd / u_tau**2
uvp_cfd = -uv_cfd / u_tau**2       # sign convention: positive in channel core

# Only plot half-channel
half = yp_cfd <= Re_tau + 5

# ── Law of the wall reference lines ──────────────────────────────────────────
yp_visc  = np.linspace(0.5, 11.5, 50)   # viscous sublayer: U+ = y+
yp_log   = np.logspace(np.log10(11.5), np.log10(Re_tau), 100)
kappa, B = 0.41, 5.2
Up_log   = (1/kappa) * np.log(yp_log) + B

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# 1. Mean velocity U+
ax = axes[0]
ax.semilogx(dns_yp,      dns_Up,      "k-",  lw=1.5, label="DNS — Moser et al. (1999)")
ax.semilogx(yp_log,      Up_log,      "k--", lw=0.9, label=f"Log law: U⁺ = (1/{kappa}) ln y⁺ + {B}")
ax.semilogx(yp_visc,     yp_visc,     "k:",  lw=0.9, label="Viscous sublayer: U⁺ = y⁺")
ax.semilogx(yp_cfd[half], Up_cfd[half], color="#e41a1c", lw=2, label="LES — WALE")
ax.set_xlabel("y⁺", fontsize=11)
ax.set_ylabel("U⁺ = ⟨U⟩ / u_τ", fontsize=11)
ax.set_title("Mean streamwise velocity", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim(0.1, Re_tau * 2)
ax.set_ylim(0, 25)

# 2. Turbulent kinetic energy k+
ax = axes[1]
ax.plot(dns_yp_s,       dns_k,       "k-",  lw=1.5, label="DNS — Moser et al. (1999)")
ax.plot(yp_cfd[half],   kp_cfd[half], color="#e41a1c", lw=2, label="LES — WALE")
ax.set_xlabel("y⁺", fontsize=11)
ax.set_ylabel("k⁺ = k / u_τ²", fontsize=11)
ax.set_title("Turbulent kinetic energy", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, Re_tau)
ax.set_ylim(0, 5)

# 3. Reynolds shear stress -<u'v'>+
ax = axes[2]
ax.plot(dns_yp_s,       -dns_uv,      "k-",  lw=1.5, label="DNS — Moser et al. (1999)")
ax.plot(yp_cfd[half],   uvp_cfd[half], color="#e41a1c", lw=2, label="LES — WALE")
ax.set_xlabel("y⁺", fontsize=11)
ax.set_ylabel("−⟨u′v′⟩⁺", fontsize=11)
ax.set_title("Reynolds shear stress", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, Re_tau)
ax.set_ylim(0, 1.2)

fig.suptitle(f"Channel flow Reτ = {Re_tau:.0f} — LES (WALE) vs DNS",
             fontsize=12, y=1.01)
fig.tight_layout()
fig.savefig(IMAGES / "validation.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("validation.png done")

# ── Convergence: bulk velocity and centerline U+ history ─────────────────────
times_all = sorted(pp_dir.iterdir(), key=lambda p: float(p.name))
t_hist, Ub_hist, Ucl_hist = [], [], []
for t_dir in times_all:
    try:
        dd = np.loadtxt(t_dir / "layerAverage.xy", comments="#")
        if dd.ndim < 2:
            continue
        t_hist.append(float(t_dir.name))
        Ub_hist.append(np.trapz(dd[:, 3], dd[:, 0]) / dd[:, 0].max())
        Ucl_hist.append(dd[:, 3].max())
    except Exception:
        continue

t_hist  = np.array(t_hist)
Ub_hist = np.array(Ub_hist)
Ucl_hist = np.array(Ucl_hist)

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t_hist, Ub_hist / u_tau,  color="#377eb8", lw=1.5, label="Bulk U⁺  (= Ū/u_τ)")
ax.plot(t_hist, Ucl_hist / u_tau, color="#e41a1c", lw=1.5, label="Centreline U⁺")
ax.axhline(dns_Up.max(), color="k", lw=1, ls="--", label=f"DNS centreline U⁺ = {dns_Up.max():.1f}")
ax.set_xlabel("Flow-through time  t·u_τ/δ  →  t / δ", fontsize=10)
ax.set_ylabel("U⁺", fontsize=10)
ax.set_title("Channel Reτ = 395 — velocity convergence history (LES / WALE)", fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(IMAGES / "convergence.png", dpi=150)
plt.close(fig)
print("convergence.png done")

# ── Velocity contour plots from VTK (foamToVTK output) ───────────────────────
VTK_FILE = Path("VTK/openfoam-channel395_5000.vtk")
if VTK_FILE.exists():
    try:
        import vtk as _vtk
        import vtk.util.numpy_support as vtk_np

        reader = _vtk.vtkUnstructuredGridReader()
        reader.SetFileName(str(VTK_FILE))
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        grid = reader.GetOutput()

        pts = vtk_np.vtk_to_numpy(grid.GetPoints().GetData())
        pd  = grid.GetPointData()

        def get_arr(name):
            a = pd.GetArray(name)
            return vtk_np.vtk_to_numpy(a) if a else None

        UMean = get_arr("UMean")          # (N,3) time-averaged velocity
        U_inst = get_arr("U")             # (N,3) instantaneous velocity

        # Extract z-midplane slice (z closest to domain centre z=1.0)
        z_unique = np.unique(pts[:, 2])
        z_mid    = z_unique[np.argmin(np.abs(z_unique - 1.0))]
        mask_z   = np.abs(pts[:, 2] - z_mid) < 1e-4

        x_sl  = pts[mask_z, 0]
        y_sl  = pts[mask_z, 1]
        Um_sl = UMean[mask_z, 0]   # streamwise component of time-averaged U

        # Wall units for y-axis label
        yp_sl = y_sl * u_tau / nu
        Up_sl = Um_sl / u_tau

        fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))

        # ── Panel 1: mean streamwise velocity UMean_x ─────────────────────────
        ax = axes[0]
        sc = ax.scatter(x_sl, y_sl, c=Up_sl, cmap="RdYlBu_r", s=4,
                        vmin=0, vmax=Up_sl.max() * 1.05, rasterized=True)
        cb = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label("U⁺ = ⟨U_x⟩ / u_τ", fontsize=9)
        ax.axhline(0, color="k", lw=1.5)
        ax.axhline(2, color="k", lw=1.5)
        ax.set_xlabel("x (m) — streamwise", fontsize=10)
        ax.set_ylabel("y (m) — wall-normal", fontsize=10)
        ax.set_title("Time-averaged streamwise velocity ⟨U_x⟩⁺\n"
                     f"z-midplane slice  (Reτ = {Re_tau:.0f}, WALE LES)", fontsize=9)
        ax.set_xlim(0, 4); ax.set_ylim(0, 2)

        # ── Panel 2: instantaneous streamwise velocity ─────────────────────────
        if U_inst is not None:
            Ui_sl  = U_inst[mask_z, 0] / u_tau
            ax2 = axes[1]
            sc2 = ax2.scatter(x_sl, y_sl, c=Ui_sl, cmap="RdYlBu_r", s=4,
                              vmin=0, vmax=Ui_sl.max() * 1.05, rasterized=True)
            cb2 = fig.colorbar(sc2, ax=ax2, pad=0.02)
            cb2.set_label("U⁺ = U_x / u_τ", fontsize=9)
            ax2.axhline(0, color="k", lw=1.5)
            ax2.axhline(2, color="k", lw=1.5)
            ax2.set_xlabel("x (m) — streamwise", fontsize=10)
            ax2.set_ylabel("y (m) — wall-normal", fontsize=10)
            ax2.set_title("Instantaneous streamwise velocity U_x⁺\n"
                          "(turbulent structures visible near walls)", fontsize=9)
            ax2.set_xlim(0, 4); ax2.set_ylim(0, 2)
        else:
            axes[1].set_visible(False)

        fig.suptitle(f"Channel Reτ = {Re_tau:.0f} — LES velocity field at t = 5000 s",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(IMAGES / "contours.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("contours.png done")

    except ImportError:
        print("vtk python package not found — skipping contours.png")
else:
    print(f"VTK file {VTK_FILE} not found — run foamToVTK first")
