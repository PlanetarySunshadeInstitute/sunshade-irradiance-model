"""
rotate
======

Rotate a sun-fixed shade field into geographic (latitude, longitude)
coordinates, given the Sun's subsolar latitude and longitude for a moment in
time.

Scientific context
------------------
The PSF irradiance model produces, for each simulated day, a 2-D dimming
factor field DF(lat, lon) in a *sun-fixed* frame: the substellar point
(directly overhead) sits at (lat=0, lon=0) on the input grid.  This matches
the MATLAB upstream's stated convention (see the top-level ``README.md``:
"The sub-stellar point is placed at CESM longitude 0°, consistent with the
noon UTC computation time used throughout the model").

To use this in a climate model we need DF in the *geographic* frame, where
the substellar point shifts ±23.44° in latitude over the year as the Sun
crosses the celestial equator.  The transformation is a rigid rotation on
the sphere: rotate the sun-fixed grid by +δ around the rotational axis that
keeps longitude 0 invariant — equivalently, "lift" the grid northward by
δ degrees while preserving its longitudinal alignment.

This module implements that rotation, expressed as a backward map: for each
*output* geographic cell (φ, θ) on the destination grid, compute the
*source* sun-fixed coordinates (φ₀, θ₀) that should be sampled, then
bilinearly interpolate the input field at those source coordinates.

The math
--------
Let

    φ      = output latitude  (geographic)
    θ      = output longitude (geographic)
    δ      = subsolar latitude  (Sun's declination, signed)
    λ      = subsolar longitude (= 0° at noon UT, the MATLAB model's
             computation time)
    φ₀, θ₀ = source latitude, longitude in the sun-fixed frame

Then the inverse rotation (backward map) is:

    sin(φ₀) = sin(φ)·cos(δ) − cos(φ)·sin(δ)·cos(θ − λ)

    y       = sin(θ − λ)·cos(φ)
    x       = cos(φ)·cos(δ)·cos(θ − λ) + sin(φ)·sin(δ)
    θ₀     = atan2(y, x)              (mod 360°)

Sanity checks built into the algorithm:

* δ = 0  → identity rotation: φ₀ ≡ φ and θ₀ ≡ θ.
* δ = 23.44° (June solstice), output cell (φ=23.44°, θ=0°) → source
  (φ₀=0, θ₀=0): the geographic substellar point pulls back to the
  sun-fixed substellar point.
* δ = −23.44° (December solstice) is the mirror image.

These are unit-tested in ``tests/test_rotate.py``.

Historical note
---------------
Earlier revisions of this module carried a ``+180°`` offset in the θ₀
formula and defaulted ``subsolar_lon_deg`` to 180°, on the assumption (from
Illeana's original Python specification) that the sun-fixed source had its
substellar point at lon=180°.  Empirically, Matt's MATLAB output places
substellar at SRC lon=0° instead — the previous assumption produced output
whose shadow landed at the antipode of the correct position.  The offset
has been removed; tests have been updated in lock-step.

Why pure numpy and not scipy
----------------------------
Illeana's reference uses ``scipy.interpolate.RegularGridInterpolator`` with
``method='linear'``.  We implement the same bilinear interpolation in plain
numpy here.  The result is mathematically identical (bilinear is bilinear),
but the algorithm is laid out explicitly in ~10 lines of code so a
scientific reader can follow exactly which weights are applied to which
neighbours.  This satisfies the project's "lots of visibility" mandate
(CLAUDE.md): no library opacity sitting between the rotation math and the
final numbers.

For longitudes near the prime meridian wrap: the source grid is padded
with a wrap-around column so that an interpolation between lon=358.75°
and lon=0° behaves correctly (the values are taken from each other) rather
than falling off the grid and being filled with the no-shade value.

Default fill value
------------------
Anywhere the backward map lands outside the source grid, the output is set
to ``fill_value`` (default 1.0 = "no shade", per Illeana).  In practice
this only happens at latitudes outside [−90, 90], which floating-point
clipping prevents — so fill should be reachable only via the optional
restricted source grid.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


logger = logging.getLogger(__name__)


# Conventions
DEFAULT_SUBSOLAR_LON_DEG = 0.0     # noon UT places substellar at longitude 0°,
                                   # matching the MATLAB upstream's stated
                                   # convention (see top-level README).
DEFAULT_FILL_VALUE = 1.0           # no shade
EARTH_OBLIQUITY_DEG = 23.44        # used only for sanity-check assertions


# -----------------------------------------------------------------------------
# Return type for the diagnostic version of rotate()
# -----------------------------------------------------------------------------

@dataclass
class RotationDiagnostics:
    """Optional per-call diagnostics returned by :func:`rotate_shade_field`
    when ``return_diagnostics=True``.  Useful for tests and for the
    transparency-style logging expected by Illeana's spec ("write an example
    of what it is taking")."""

    nlats: int
    nlons: int
    subsolar_lat_deg: float
    subsolar_lon_deg: float
    n_filled: int            # cells whose backward map landed outside grid
    input_min: float
    input_max: float
    output_min: float
    output_max: float


# -----------------------------------------------------------------------------
# Default source grid — matches Illeana's spec
# -----------------------------------------------------------------------------

def default_source_grid(nlats: int, nlons: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    The lat/lon grid Illeana's reference code assumes:

        lats = np.linspace(-90, 90, nlats)        # poles included
        lons = np.linspace(0, 360, nlons)         # 0 and 360 both included

    Returns ``(lats, lons)``, both 1-D float64.

    Note
    ----
    The CESM "f09" grid (1° × 1.25°) actually uses
    ``lons = np.linspace(0, 360, 288, endpoint=False)`` — i.e. 0° to
    358.75° at exact 1.25° spacing — and the pole rows are at ±90° on a
    192-point lat grid.  When that is the right convention for the input
    file, pass it explicitly via ``source_lats_deg`` / ``source_lons_deg``.
    """
    return (
        np.linspace(-90.0, 90.0, nlats),
        np.linspace(0.0, 360.0, nlons),
    )


# -----------------------------------------------------------------------------
# Bilinear interpolation on a regular grid (lat, lon)
# -----------------------------------------------------------------------------

def _bilinear_interp_with_lon_wrap(
    field: np.ndarray,
    src_lats: np.ndarray,
    src_lons: np.ndarray,
    query_lats: np.ndarray,
    query_lons: np.ndarray,
    fill_value: float,
) -> np.ndarray:
    """
    Bilinear interpolation of ``field`` at the (query_lat, query_lon) points,
    with longitude treated periodically.

    Parameters
    ----------
    field        : (nlat, nlon) source field, indexed [lat_idx, lon_idx].
    src_lats     : (nlat,) strictly increasing latitude axis (degrees).
    src_lons     : (nlon,) strictly increasing longitude axis (degrees).
                   May be 0..360 with both endpoints (Illeana's default), or
                   0..360 with right endpoint excluded (CESM f09 style).
    query_lats   : array of query latitudes (degrees), any shape.
    query_lons   : array of query longitudes (degrees), same shape as
                   ``query_lats``.  Will be normalised modulo 360.
    fill_value   : value to write where the query latitude is outside the
                   source latitude range.  Longitude is periodic so longitude
                   never falls "outside".

    Returns
    -------
    out : ndarray shaped like ``query_lats`` with the interpolated values.

    Implementation
    --------------
    Latitude: standard linear interpolation between the two enclosing rows,
    with weights set by fractional distance to the lower row.

    Longitude: extend the source axis by appending the column ``lon =
    src_lons[0] + 360`` (with values copied from column 0) so that an
    interpolation between ``src_lons[-1]`` and the wrap-around ``360 + 0`` is
    well-defined.  This avoids artefacts in a small wedge at the dateline /
    prime meridian when ``src_lons`` does not include both 0 and 360.

    All operations are vectorised over the query points.
    """
    nlat = src_lats.size
    nlon = src_lons.size
    if field.shape != (nlat, nlon):
        raise ValueError(
            f"field shape {field.shape} does not match grid "
            f"({nlat} lat × {nlon} lon)"
        )

    # --- Pad the source grid with one wrap-around longitude column ---------
    # New lon axis: src_lons concatenated with (src_lons[0] + 360).
    # New field:    field with column 0 appended on the right.
    # This is cheap (one extra column) and lets us treat lon as a closed
    # interval [src_lons[0], src_lons[0] + 360] for interpolation, while
    # remaining mathematically identical to a periodic interpolation.
    lons_ext = np.concatenate([src_lons, [src_lons[0] + 360.0]])
    field_ext = np.concatenate([field, field[:, :1]], axis=1)
    nlon_ext = lons_ext.size  # = nlon + 1

    # --- Normalise query longitudes into [src_lons[0], src_lons[0] + 360) --
    q_lon_norm = ((query_lons - src_lons[0]) % 360.0) + src_lons[0]

    # --- Find the bracketing latitude indices ------------------------------
    # np.searchsorted gives the insertion index that would keep src_lats
    # sorted; we then clamp to a valid bracket [i_lo, i_hi=i_lo+1].
    i_hi = np.searchsorted(src_lats, query_lats, side="right")
    in_lat_range = (query_lats >= src_lats[0]) & (query_lats <= src_lats[-1])
    i_hi = np.clip(i_hi, 1, nlat - 1)
    i_lo = i_hi - 1

    lat_lo = src_lats[i_lo]
    lat_hi = src_lats[i_hi]
    # Avoid division by zero in the (impossible for monotonic grids) edge
    # where lat_lo == lat_hi.
    denom_lat = np.where(lat_hi != lat_lo, lat_hi - lat_lo, 1.0)
    w_lat = (query_lats - lat_lo) / denom_lat

    # --- Find the bracketing longitude indices on the *extended* axis ------
    j_hi = np.searchsorted(lons_ext, q_lon_norm, side="right")
    j_hi = np.clip(j_hi, 1, nlon_ext - 1)
    j_lo = j_hi - 1
    lon_lo = lons_ext[j_lo]
    lon_hi = lons_ext[j_hi]
    denom_lon = np.where(lon_hi != lon_lo, lon_hi - lon_lo, 1.0)
    w_lon = (q_lon_norm - lon_lo) / denom_lon

    # --- Gather the four corner values and combine -------------------------
    f00 = field_ext[i_lo, j_lo]   # (lat_lo, lon_lo)
    f01 = field_ext[i_lo, j_hi]   # (lat_lo, lon_hi)
    f10 = field_ext[i_hi, j_lo]   # (lat_hi, lon_lo)
    f11 = field_ext[i_hi, j_hi]   # (lat_hi, lon_hi)

    out = (
        (1.0 - w_lat) * (1.0 - w_lon) * f00
        + (1.0 - w_lat) * w_lon       * f01
        + w_lat        * (1.0 - w_lon) * f10
        + w_lat        * w_lon        * f11
    )

    # --- Apply fill where query latitude is outside the source range -------
    out = np.where(in_lat_range, out, fill_value)

    return out


# -----------------------------------------------------------------------------
# Backward map (the inverse rotation)
# -----------------------------------------------------------------------------

def backward_map(
    out_lat_deg: np.ndarray,
    out_lon_deg: np.ndarray,
    subsolar_lat_deg: float,
    subsolar_lon_deg: float = DEFAULT_SUBSOLAR_LON_DEG,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    For each *output* geographic cell (φ, θ), return the *source* sun-fixed
    coordinates (φ₀, θ₀) to sample.

    Vectorised: ``out_lat_deg`` and ``out_lon_deg`` may be any shape (so long
    as they broadcast); the returned arrays have that broadcast shape.

    Returns
    -------
    src_lat_deg, src_lon_deg : two ndarrays in degrees.  ``src_lon_deg`` is
        normalised to ``[0, 360)``.
    """
    phi    = np.radians(out_lat_deg)
    theta  = np.radians(out_lon_deg)
    delta  = np.radians(subsolar_lat_deg)
    lambd  = np.radians(subsolar_lon_deg)

    # sin(φ₀) — clip to avoid arcsin domain errors from float roundoff at poles
    sin_phi0 = (
        np.sin(phi) * np.cos(delta)
        - np.cos(phi) * np.sin(delta) * np.cos(theta - lambd)
    )
    phi0 = np.arcsin(np.clip(sin_phi0, -1.0, 1.0))

    # θ₀ via atan2.  No additive offset: the sun-fixed grid places its
    # substellar point at lon=0°, matching the MATLAB upstream, so the
    # geographic substellar (δ, λ) pulls back to source (0°, 0°) cleanly.
    y = np.sin(theta - lambd) * np.cos(phi)
    x = (
        np.cos(phi) * np.cos(delta) * np.cos(theta - lambd)
        + np.sin(phi) * np.sin(delta)
    )
    theta0 = np.arctan2(y, x)

    src_lat_deg = np.degrees(phi0)
    src_lon_deg = np.degrees(theta0) % 360.0

    return src_lat_deg, src_lon_deg


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

def rotate_shade_field(
    shade_data: np.ndarray,
    subsolar_lat_deg: float,
    subsolar_lon_deg: float = DEFAULT_SUBSOLAR_LON_DEG,
    source_lats_deg: Optional[np.ndarray] = None,
    source_lons_deg: Optional[np.ndarray] = None,
    fill_value: float = DEFAULT_FILL_VALUE,
    return_diagnostics: bool = False,
):
    """
    Rotate a sun-fixed shade field into geographic coordinates.

    Parameters
    ----------
    shade_data : 2-D ndarray of shape ``(nlats, nlons)``
        Dimming factor in the sun-fixed frame.  Axis 0 = latitude
        (south-to-north), axis 1 = longitude (0° eastward to 360°).  The
        substellar point sits at (lat=0, lon=0) on this grid, matching the
        MATLAB upstream's convention.
    subsolar_lat_deg : float
        Sun's declination δ for this snapshot, in signed decimal degrees
        (positive = northern hemisphere).
    subsolar_lon_deg : float
        Sun's geographic longitude at the snapshot moment.  Defaults to 0°
        (noon UT, when the substellar point is over the prime meridian —
        the convention recorded in the top-level ``README.md``).
    source_lats_deg, source_lons_deg : optional 1-D ndarrays
        The latitude / longitude axis of ``shade_data``.  If omitted,
        Illeana's defaults are used: ``linspace(-90, 90, nlats)`` and
        ``linspace(0, 360, nlons)``.
    fill_value : float, default 1.0
        Value written for any output cell whose backward map lands outside
        the source grid.  ``1.0`` = "no shade" — the conservative choice for
        an irradiance multiplier.
    return_diagnostics : bool, default False
        If True, return ``(rotated, RotationDiagnostics(...))`` instead of
        just ``rotated``.  Useful for transparency logging.

    Returns
    -------
    rotated : ndarray of the same shape as ``shade_data``, in geographic
        coordinates on the same lat/lon axes.
    diagnostics : RotationDiagnostics (only if ``return_diagnostics=True``)

    Raises
    ------
    ValueError
        If ``shade_data`` is not 2-D, if the supplied source grid does not
        match its shape, or if ``subsolar_lat_deg`` is outside [−90, 90].
    """
    # --- Validation --------------------------------------------------------
    arr = np.asarray(shade_data)
    if arr.ndim != 2:
        raise ValueError(
            f"shade_data must be 2-D (lat, lon); got shape {arr.shape}"
        )
    nlats, nlons = arr.shape

    if not np.isfinite(subsolar_lat_deg):
        raise ValueError(f"subsolar_lat_deg must be finite; got {subsolar_lat_deg!r}")
    if not (-90.0 <= subsolar_lat_deg <= 90.0):
        raise ValueError(
            f"subsolar_lat_deg must lie in [-90, 90]; got {subsolar_lat_deg}"
        )
    # An alarm bell: physically, the Sun's declination is bounded by Earth's
    # obliquity (~23.44°). Anything more than 25° is almost certainly an
    # input mistake (e.g. radians passed as degrees).  Warn but don't fail.
    if abs(subsolar_lat_deg) > EARTH_OBLIQUITY_DEG + 1.5:
        logger.warning(
            "subsolar_lat_deg = %+.3f° exceeds Earth's obliquity (±%.2f°); "
            "input may be wrong (radians vs. degrees?).",
            subsolar_lat_deg, EARTH_OBLIQUITY_DEG,
        )

    if source_lats_deg is None:
        src_lats, src_lons_default = default_source_grid(nlats, nlons)
    else:
        src_lats = np.asarray(source_lats_deg, dtype=np.float64)
    if source_lons_deg is None:
        if source_lats_deg is None:
            src_lons = src_lons_default
        else:
            _, src_lons = default_source_grid(nlats, nlons)
    else:
        src_lons = np.asarray(source_lons_deg, dtype=np.float64)

    if src_lats.shape != (nlats,):
        raise ValueError(
            f"source_lats_deg shape {src_lats.shape} does not match "
            f"shade_data lat axis ({nlats},)"
        )
    if src_lons.shape != (nlons,):
        raise ValueError(
            f"source_lons_deg shape {src_lons.shape} does not match "
            f"shade_data lon axis ({nlons},)"
        )
    if not np.all(np.diff(src_lats) > 0):
        raise ValueError("source_lats_deg must be strictly increasing")
    if not np.all(np.diff(src_lons) > 0):
        raise ValueError("source_lons_deg must be strictly increasing")

    # --- Build the output grid (same axes as input) -----------------------
    # The output cells are the same lat/lon points as the input grid.  We
    # use meshgrid with indexing='xy' so lat_grid[i, j] varies with i (lat)
    # and lon_grid[i, j] varies with j (lon) — consistent with the (lat,
    # lon) shape of the field.
    lon_grid, lat_grid = np.meshgrid(src_lons, src_lats)   # both (nlats, nlons)

    # --- Backward map every output cell to its source location ------------
    src_lat_q, src_lon_q = backward_map(
        lat_grid, lon_grid, subsolar_lat_deg, subsolar_lon_deg,
    )

    # --- Bilinearly interpolate the input field at those source points ----
    rotated = _bilinear_interp_with_lon_wrap(
        field=arr,
        src_lats=src_lats,
        src_lons=src_lons,
        query_lats=src_lat_q,
        query_lons=src_lon_q,
        fill_value=fill_value,
    )

    if return_diagnostics:
        # Count cells filled from the fill value (i.e. backward map landed
        # outside the source latitude range).  In practice this should be 0
        # for a global grid covering [-90, 90].
        in_range = (src_lat_q >= src_lats[0]) & (src_lat_q <= src_lats[-1])
        diag = RotationDiagnostics(
            nlats=nlats, nlons=nlons,
            subsolar_lat_deg=float(subsolar_lat_deg),
            subsolar_lon_deg=float(subsolar_lon_deg),
            n_filled=int(np.sum(~in_range)),
            input_min=float(arr.min()),
            input_max=float(arr.max()),
            output_min=float(rotated.min()),
            output_max=float(rotated.max()),
        )
        return rotated, diag

    return rotated
