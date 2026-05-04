"""
plot_annual_mean_df.py
======================
Figure 5 for the PSI methods writeup.

Reads the 31-year geographic-frame dimming-factor NetCDF
(L1-RX-2026-04-17-001_df_2035_2065.nc), computes the time-mean DF over
all 11,323 daily snapshots, and saves a Mercator (PlateCarree) filled
colour map with coastlines to writeup/figures/fig_mercator_df.png.

Usage
-----
    cd /path/to/repo/root
    python python/coordinate_transformation/plot_annual_mean_df.py

Dependencies
------------
    numpy, netCDF4, matplotlib, cartopy
    Install: pip install numpy netCDF4 matplotlib cartopy

    If cartopy is unavailable, the script falls back to a plain
    matplotlib pcolormesh (no coastlines) and prints a warning.
"""

from __future__ import annotations

import os
import sys
import warnings
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")           # non-interactive backend — safe for scripts
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Paths — all relative to the repo root (the directory above python/)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT  = os.path.abspath(os.path.join(_SCRIPT_DIR, "..", ".."))

NC_PATH   = os.path.join(_REPO_ROOT,
                         "numerical control", "exports",
                         "L1-RX-2026-04-17-001_df_2035_2065.nc")
SAVE_PATH = os.path.join(_REPO_ROOT, "writeup", "figures", "fig_mercator_df.png")


# ---------------------------------------------------------------------------
# Read NC file
# ---------------------------------------------------------------------------
print(f"Reading NC file:\n  {NC_PATH}")

if not os.path.isfile(NC_PATH):
    sys.exit(f"ERROR: NC file not found:\n  {NC_PATH}\n"
             "Adjust NC_PATH in this script if the file has moved.")

with nc.Dataset(NC_PATH) as ds:
    lat = ds["latitude"][:]          # (192,)  degrees_north, south-to-north
    lon = ds["longitude"][:]         # (288,)  degrees_east,  0 – 358.75°
    df  = ds["DF"][:]                # (11323, 192, 288)  masked array

print(f"  DF shape : {df.shape}")
print(f"  lat      : {float(lat[0]):.3f} … {float(lat[-1]):.3f}°")
print(f"  lon      : {float(lon[0]):.3f} … {float(lon[-1]):.3f}°")


# ---------------------------------------------------------------------------
# Time-mean
# ---------------------------------------------------------------------------
print("Computing time-mean DF …")
df_mean = np.ma.mean(df, axis=0)    # (192, 288)
print(f"  mean range: {float(df_mean.min()):.6f} – {float(df_mean.max()):.6f}")


# ---------------------------------------------------------------------------
# Colormap — blue (deep shading) → white (no shading)
# Consistent with Figures 3 and 4 (same RGB interpolation as the MATLAB
# irradiance diagnostic: R=[0.10,1.00], G=[0.30,1.00], B=[0.80,1.00]).
# ---------------------------------------------------------------------------
n_colors = 256
frac     = np.linspace(0.0, 1.0, n_colors)
r_ch     = np.interp(frac, [0.0, 1.0], [0.10, 1.00])
g_ch     = np.interp(frac, [0.0, 1.0], [0.30, 1.00])
b_ch     = np.interp(frac, [0.0, 1.0], [0.80, 1.00])
cmap_psi = mcolors.ListedColormap(np.column_stack([r_ch, g_ch, b_ch]))


# ---------------------------------------------------------------------------
# Close the CESM longitude seam (0° … 358.75° → add column at 360°)
# ---------------------------------------------------------------------------
try:
    from cartopy.util import add_cyclic_point
    df_plot, lon_plot = add_cyclic_point(df_mean, coord=lon)
    _has_cyclic = True
except Exception:
    df_plot  = df_mean
    lon_plot = lon
    _has_cyclic = False


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    _has_cartopy = True
except ImportError:
    _has_cartopy = False
    warnings.warn(
        "cartopy is not installed — falling back to plain matplotlib "
        "(no coastlines).  Install with: pip install cartopy",
        stacklevel=1,
    )

fig = plt.figure(figsize=(11, 5.5), facecolor="white")

if _has_cartopy:
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.PlateCarree(central_longitude=180))
    pcm = ax.pcolormesh(
        lon_plot, lat, df_plot,
        cmap=cmap_psi,
        vmin=float(df_mean.min()),
        vmax=1.0,
        transform=ccrs.PlateCarree(),
        rasterized=True,
    )
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6, edgecolor="black", zorder=3)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.3, edgecolor="0.4",   zorder=3)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.35,
        color="gray",
        linestyle="--",
        alpha=0.7,
    )
    gl.top_labels   = False
    gl.right_labels = False
    ax.set_global()
else:
    # Plain fallback — no coastlines
    ax = fig.add_subplot(1, 1, 1)
    pcm = ax.pcolormesh(lon_plot, lat, df_plot, cmap=cmap_psi,
                        vmin=float(df_mean.min()), vmax=1.0, rasterized=True)
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°E)", fontsize=10)
    ax.set_ylabel("Latitude (°N)", fontsize=10)

# Colorbar
cb = fig.colorbar(pcm, ax=ax, orientation="vertical",
                  pad=0.03, fraction=0.025, shrink=0.85)
cb.set_label("Annual-mean dimming factor", fontsize=10)
cb.ax.tick_params(labelsize=9)

ax.set_title(
    "Annual-mean dimming factor, geographic frame  (2035–2065)",
    fontsize=12, fontweight="bold", pad=10,
)

fig.tight_layout()


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(SAVE_PATH), exist_ok=True)
fig.savefig(SAVE_PATH, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved figure:\n  {SAVE_PATH}")
