"""
plot_rotation_comparison.py
===========================
Figure 7 for the PSI methods writeup — four seasonal snapshots.

Generates one figure per cardinal date (March equinox, June solstice,
September equinox, December solstice), each with two side-by-side panels:

  Left  — Dimming factor in the sun-fixed Earth-disc frame.
           Plotted in (sin θ, sin φ) space so the disc appears circular.
           Sub-stellar point always centred at (0°, 0°).

  Right — The same day's DF in the geographic frame after the Python
           coordinate transformation.  Sub-stellar point at (0°E, δ°N/S)
           where δ is the solar declination for that date.

All four figures share the same colormap (PSI blue → white) and the same
colour scale (vmin = global minimum across all four days).

Output files
------------
    writeup/figures/fig_rotation_comparison_mar_equinox.png
    writeup/figures/fig_rotation_comparison_jun_solstice.png
    writeup/figures/fig_rotation_comparison_sep_equinox.png
    writeup/figures/fig_rotation_comparison_dec_solstice.png

Usage
-----
    cd /path/to/repo/root
    python python/coordinate_transformation/plot_rotation_comparison.py

Dependencies
------------
    numpy, netCDF4, matplotlib, cartopy
    Install: pip install numpy netCDF4 matplotlib cartopy
    cartopy is required for the right panel.  If unavailable the script
    falls back to a plain matplotlib axes (no coastlines) with a warning.
"""

from __future__ import annotations

import os
import sys
import warnings
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT  = os.path.abspath(os.path.join(_SCRIPT_DIR, "..", ".."))

NC_DISC   = os.path.join(_REPO_ROOT, "numerical control", "exports",
                         "L1-RX-2026-04-17-001.nc")
NC_GEO    = os.path.join(_REPO_ROOT, "numerical control", "exports",
                         "L1-RX-2026-04-17-001_df_2035_2065.nc")
SAVE_DIR  = os.path.join(_REPO_ROOT, "writeup", "figures")

for p in (NC_DISC, NC_GEO):
    if not os.path.isfile(p):
        sys.exit(f"ERROR: file not found:\n  {p}")


# ---------------------------------------------------------------------------
# Cardinal dates
#   doy        : day-of-year (1-indexed) in the 365-day sun-fixed file
#   idx_geo    : 0-indexed row in the geographic file (starts 2035-01-01)
#   decl_deg   : approximate solar declination on that date (degrees)
#   date_str   : human-readable date label
#   short_name : used in the output filename
# ---------------------------------------------------------------------------
DATES = [
    dict(doy=80,  idx_geo=79,  decl_deg=  0.0,
         date_str="2035-03-21", label="March equinox",
         short_name="mar_equinox"),
    dict(doy=172, idx_geo=171, decl_deg=+23.44,
         date_str="2035-06-21", label="June solstice",
         short_name="jun_solstice"),
    dict(doy=266, idx_geo=265, decl_deg=  0.0,
         date_str="2035-09-23", label="September equinox",
         short_name="sep_equinox"),
    dict(doy=355, idx_geo=354, decl_deg=-23.44,
         date_str="2035-12-21", label="December solstice",
         short_name="dec_solstice"),
]

SUBSOLAR_LON = 0.0   # always 0° — CESM noon-UTC convention


# ---------------------------------------------------------------------------
# Read sun-fixed disc file (full, all 365 days)
#   MATLAB stores (lon=288, lat=192, time=365) in Fortran order;
#   Python netCDF4 reads it as (time=365, lat=192, lon=288).
#
#   Disc reconstruction (view_NC_File___0___0.m):
#     CESM lon cols 217–288 (0-indexed 216–287) → theta −89.375° to −0.625° (west)
#     CESM lon cols   1– 72 (0-indexed   0– 71) → theta  +0.625° to +89.375° (east)
#   Combined: 144 columns, theta = linspace(−89.375, 89.375, 144)
# ---------------------------------------------------------------------------
print(f"Reading sun-fixed disc file …\n  {NC_DISC}")
disc_cols = list(range(216, 288)) + list(range(0, 72))   # 144 CESM lon indices
theta_deg = np.linspace(-89.375, 89.375, 144)

with nc.Dataset(NC_DISC) as ds:
    phi_deg    = ds["lat"][:]          # (192,)
    df_disc_4  = np.stack(
        [ds["df"][d["doy"] - 1, :, :][:, disc_cols] for d in DATES],
        axis=0,
    )                                  # (4, 192, 144)  phi × theta, one slab per date

print(f"  Loaded 4 disc slabs, shape: {df_disc_4.shape}")


# ---------------------------------------------------------------------------
# Read geographic file (only the 4 days we need)
# ---------------------------------------------------------------------------
print(f"Reading geographic frame file …\n  {NC_GEO}")

with nc.Dataset(NC_GEO) as ds:
    lat_geo   = ds["latitude"][:]      # (192,)
    lon_geo   = ds["longitude"][:]     # (288,)  0 – 358.75°
    df_geo_4  = np.stack(
        [ds["DF"][d["idx_geo"], :, :] for d in DATES],
        axis=0,
    )                                  # (4, 192, 288)

print(f"  Loaded 4 geo slabs, shape: {df_geo_4.shape}")


# ---------------------------------------------------------------------------
# Shared colour scale — consistent vmin across all four dates and both panels
# ---------------------------------------------------------------------------
vmin = min(float(df_disc_4.min()), float(np.ma.min(df_geo_4)))
vmax = 1.0
print(f"  Shared colour scale: [{vmin:.6f}, {vmax:.6f}]")


# ---------------------------------------------------------------------------
# PSI colormap: blue (deep shading) → white (no shading)
# Same RGB ramp as MATLAB diagnostic figures (Figs 3 and 4).
# ---------------------------------------------------------------------------
n_colors = 256
frac     = np.linspace(0.0, 1.0, n_colors)
cmap_psi = mcolors.ListedColormap(
    np.column_stack([
        np.interp(frac, [0, 1], [0.10, 1.00]),
        np.interp(frac, [0, 1], [0.30, 1.00]),
        np.interp(frac, [0, 1], [0.80, 1.00]),
    ])
)


# ---------------------------------------------------------------------------
# Cartopy — try once; fall back gracefully
# ---------------------------------------------------------------------------
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.util import add_cyclic_point
    _has_cartopy = True
    print("  cartopy available — coastlines will be drawn.")
except ImportError:
    _has_cartopy = False
    warnings.warn(
        "cartopy not installed — right panels will have no coastlines.\n"
        "Install with: pip install cartopy",
        stacklevel=1,
    )


# ---------------------------------------------------------------------------
# Reusable geometry for the disc panel
# ---------------------------------------------------------------------------
t_circ    = np.linspace(0, 2 * np.pi, 500)
tick_deg  = [-90, -60, -30, 0, 30, 60, 90]
tick_sin  = np.sin(np.deg2rad(tick_deg))
tick_lbl  = [f"{d}°" for d in tick_deg]
sin_theta = np.sin(np.deg2rad(theta_deg))   # (144,)
sin_phi   = np.sin(np.deg2rad(phi_deg))     # (192,)
grid_degs = [-60, -30, 0, 30, 60]

os.makedirs(SAVE_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Loop — one figure per date
# ---------------------------------------------------------------------------
for i, d in enumerate(DATES):

    disc_data = df_disc_4[i]          # (192, 144)
    geo_data  = df_geo_4[i]           # (192, 288)
    decl      = d["decl_deg"]
    decl_str  = (f"δ ≈ {decl:+.1f}°" if abs(decl) > 0.5
                 else "δ ≈ 0°  (equinox)")
    save_path = os.path.join(SAVE_DIR,
                             f"fig_rotation_comparison_{d['short_name']}.png")

    print(f"\n[{i+1}/4] {d['label']}  ({d['date_str']},  {decl_str})")

    fig = plt.figure(figsize=(13, 5.5), facecolor="white")

    # ======================================================================
    # LEFT PANEL — sun-fixed disc frame
    # ======================================================================
    ax_disc = fig.add_subplot(1, 2, 1)

    # Grey exterior, white disc interior
    ax_disc.set_facecolor([0.93, 0.93, 0.93])
    circ_bg = plt.Polygon(list(zip(np.cos(t_circ), np.sin(t_circ))),
                          facecolor="white", edgecolor="none", zorder=0)
    ax_disc.add_patch(circ_bg)

    # Data — imshow in (sin θ, sin φ) space; clipped to unit circle
    im = ax_disc.imshow(
        disc_data,
        origin="lower",
        extent=[sin_theta[0], sin_theta[-1], sin_phi[0], sin_phi[-1]],
        aspect="equal",
        cmap=cmap_psi,
        vmin=vmin, vmax=vmax,
        interpolation="nearest",
        zorder=1,
    )
    im.set_clip_path(mpatches.Circle((0, 0), 1.0,
                                     transform=ax_disc.transData))

    # Grid lines (clipped chords)
    for gd in grid_degs:
        s    = np.sin(np.deg2rad(gd))
        half = np.sqrt(max(1 - s**2, 0))
        ax_disc.plot([-half, half], [s, s],
                     color=[0.5, 0.5, 0.5], lw=0.5, ls="--", zorder=3)
        ax_disc.plot([s, s], [-half, half],
                     color=[0.5, 0.5, 0.5], lw=0.5, ls="--", zorder=3)

    # Disc boundary
    ax_disc.plot(np.cos(t_circ), np.sin(t_circ),
                 "k-", linewidth=1.5, zorder=5)

    # Sub-stellar point (always at centre in disc frame)
    ax_disc.plot(0, 0, "ko", markersize=6, zorder=6)
    ax_disc.annotate("Sub-stellar\n(0°, 0°)",
                     xy=(0, 0), xytext=(0.08, 0.10),
                     fontsize=8,
                     arrowprops=dict(arrowstyle="-", color="black", lw=0.7))

    ax_disc.set_xlim(-1.08, 1.08);  ax_disc.set_ylim(-1.08, 1.08)
    ax_disc.set_xticks(tick_sin);   ax_disc.set_xticklabels(tick_lbl, fontsize=8)
    ax_disc.set_yticks(tick_sin);   ax_disc.set_yticklabels(tick_lbl, fontsize=8)
    ax_disc.set_xlabel(r"$\theta$ — degrees from sub-stellar point", fontsize=10)
    ax_disc.set_ylabel(r"$\phi$ — degrees from sub-stellar point",   fontsize=10)
    ax_disc.set_title(f"Sun-fixed disc frame\n(day {d['doy']}, {d['label']})",
                      fontsize=11, fontweight="bold")
    ax_disc.set_aspect("equal")
    ax_disc.tick_params(direction="out")

    # ======================================================================
    # RIGHT PANEL — geographic frame
    # ======================================================================
    if _has_cartopy:
        ax_geo = fig.add_subplot(
            1, 2, 2, projection=ccrs.PlateCarree(central_longitude=0))
        geo_cyc, lon_cyc = add_cyclic_point(geo_data, coord=lon_geo)
        pcm = ax_geo.pcolormesh(
            lon_cyc, lat_geo, geo_cyc,
            cmap=cmap_psi, vmin=vmin, vmax=vmax,
            transform=ccrs.PlateCarree(), rasterized=True,
        )
        ax_geo.add_feature(cfeature.COASTLINE, linewidth=0.6,
                           edgecolor="black", zorder=4)
        gl = ax_geo.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=0.35, color="gray",
                              linestyle="--", alpha=0.7)
        gl.top_labels   = False
        gl.right_labels = False
        ax_geo.set_global()

        # Sub-stellar point
        ax_geo.plot(SUBSOLAR_LON, decl, "ko", markersize=6,
                    transform=ccrs.PlateCarree(), zorder=5)
        lon_label = f"{SUBSOLAR_LON:.0f}°E"
        lat_label = f"{abs(decl):.1f}°{'N' if decl >= 0 else 'S'}" if abs(decl) > 0.5 else "0°"
        ax_geo.annotate(
            f"Sub-stellar\n({lon_label}, {lat_label})",
            xy=ax_geo.projection.transform_point(
                SUBSOLAR_LON, decl, ccrs.PlateCarree()),
            xytext=(25, 15), textcoords="offset points",
            fontsize=8,
            arrowprops=dict(arrowstyle="-", color="black", lw=0.7),
            zorder=6,
        )
    else:
        ax_geo = fig.add_subplot(1, 2, 2)
        pcm = ax_geo.pcolormesh(lon_geo, lat_geo, geo_data,
                                cmap=cmap_psi, vmin=vmin, vmax=vmax,
                                rasterized=True)
        ax_geo.plot(SUBSOLAR_LON, decl, "ko", markersize=6)
        ax_geo.set_xlabel("Longitude (°E)", fontsize=10)
        ax_geo.set_ylabel("Latitude (°N)", fontsize=10)
        ax_geo.set_xlim(0, 360);  ax_geo.set_ylim(-90, 90)

    ax_geo.set_title(
        f"Geographic frame\n({d['date_str']},  {decl_str})",
        fontsize=11, fontweight="bold",
    )

    # ======================================================================
    # Shared colorbar
    # ======================================================================
    cbar_ax = fig.add_axes([0.92, 0.12, 0.018, 0.74])
    cb = fig.colorbar(pcm, cax=cbar_ax)
    cb.set_label("Dimming factor", fontsize=10, labelpad=8)
    cb.ax.tick_params(labelsize=9)

    fig.subplots_adjust(left=0.07, right=0.90, wspace=0.28)

    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {save_path}")

print("\nDone — 4 figures written to writeup/figures/")
