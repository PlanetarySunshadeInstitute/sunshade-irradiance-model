"""
assemble
========

Top-level orchestrator: reads Matt's 2-year sun-fixed dimming-factor NetCDF,
reads the 31-year JPL Horizons ephemeris, rotates each day's sun-fixed slab
into geographic coordinates using the Sun's declination, and streams the
result into a 30-year-plus DF(time, lat, lon) NetCDF ready for downstream
multiplication against the CMIP6 SolarForcing TSI file.

Pipeline
--------

    +-------------------------+       +--------------------------+
    | Matt's source NetCDF    |       | JPL Horizons text file   |
    |  df  (365, lat, lon)    |       |  11,323 daily (date, δ)  |
    |  dfl (366, lat, lon)    |       |  2035-01-01 .. 2065-12-31|
    |  lat, lon               |       +----------+---------------+
    +-----------+-------------+                  |
                |                                |
                v                                v
       SunFixedSourceProvider             HorizonsEphemeris
                |                                |
                +---------+----------------------+
                          |
                          v
           +---------------------------------+
           | assemble_into_writer(...)       |
           |   for t in 0..nt-1:             |
           |     is_leap = year(t) is leap   |
           |     doy     = day-of-year(t)    |
           |     slab    = dfl[doy] if leap  |
           |               else df[doy]      |
           |     rot     = rotate_shade_field|
           |               (slab, δ(t))      |
           |     writer.write_day(t, rot)    |
           +-----------------+---------------+
                             |
                             v
                +-----------------------------+
                |  Output NetCDF              |
                |    DF(time=11323, lat=192,  |
                |       lon=288)              |
                |    date(time)               |
                |    time(time)               |
                |    lat, lon                 |
                +-----------------------------+

Design notes
------------
The loop itself is pure Python + numpy — no NetCDF library — so that it can
be tested with mock providers and writers without needing netCDF4/h5py
installed.  The actual file I/O lives in two thin wrappers at the bottom of
the module:

* :class:`NetCDFSourceReader` reads Matt's source file (df, dfl, lat, lon)
  using ``netCDF4``.  Imported lazily inside ``__init__`` so the rest of the
  module is usable without netCDF4 present.
* The writer is a separate module (``write_cesm_nc``) that exposes a class
  matching :class:`WriterProtocol` below.  This file calls it only through
  that protocol, so substituting a different writer (or a mock, for tests)
  is a one-liner.

Transparency
------------
A scientist reading ``assemble_into_writer`` should be able to see, in order:

1. Every input's shape, dtype, and meaning (logged at INFO).
2. The day-zero mapping: "Horizons date D → source index I → leap? Y/N" —
   this pinpoints the one subtle assumption in the code (that ``df[0]`` is
   the Jan-1 shade pattern).  Flip that assumption and the whole pipeline
   shifts by one day.
3. Per-day progress every N days so long runs don't look hung.
4. A summary at the end: runtime, rotation min/max/mean, filled-cell count.

If any of the per-day rotations produces a physically implausible value
(e.g. DF > 1 + ε, DF < 0), the run fails loudly rather than continuing to
produce garbage.
"""

from __future__ import annotations

import calendar
import logging
import time as _time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Protocol, Tuple, Union

import numpy as np

from . import horizons_parser
from .horizons_parser import HorizonsEphemeris
from .rotate import rotate_shade_field


logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Protocols for the two I/O seams — source provider and writer.
# Both can be swapped for mocks in tests.
# -----------------------------------------------------------------------------

class SourceProviderProtocol(Protocol):
    """Anything that can yield a sun-fixed 2D shade slab for a given date."""

    source_lats_deg: np.ndarray
    source_lons_deg: np.ndarray

    def slice_for_date(self, dt: datetime) -> np.ndarray:
        """Return the sun-fixed shade slab for date ``dt``.

        Must return a 2-D float array of shape (nlat, nlon) matching
        ``source_lats_deg`` and ``source_lons_deg``.  The shape is the
        same for every date (leap or standard year).
        """
        ...


class WriterProtocol(Protocol):
    """Anything that can accept daily rotated slabs and persist them."""

    def write_day(self, t_idx: int, df_slab: np.ndarray) -> None:
        """Write the rotated slab for output time index ``t_idx`` (0-based).

        ``df_slab`` is 2-D, shape (nlat, nlon), float64.
        """
        ...

    def close(self) -> None:
        """Flush and close any underlying resources."""
        ...


# -----------------------------------------------------------------------------
# Sun-fixed source provider — implements the leap-year / day-of-year lookup.
# -----------------------------------------------------------------------------

@dataclass
class SunFixedSourceProvider:
    """
    Holds the sun-fixed shade fields for a "standard" year (365 days) and a
    "leap" year (366 days), plus the lat/lon grid.  Looks up the right slab
    for a given ``datetime`` by computing day-of-year and checking leapness.

    Convention
    ----------
    ``df_standard[0]`` and ``df_leap[0]`` are both interpreted as the
    shade pattern for **January 1** (day-of-year 1 in the standard calendar
    sense, 0-indexed here).  ``df_standard[364]`` is Dec 31.
    ``df_leap[365]`` is Dec 31 of a leap year.

    If Matt's MATLAB code has a different ordering convention (e.g.
    day-of-year starts from some other reference day), this is the single
    place to adjust it.
    """

    df_standard: np.ndarray             # (365, nlat, nlon)
    df_leap: np.ndarray                 # (366, nlat, nlon)
    source_lats_deg: np.ndarray         # (nlat,)
    source_lons_deg: np.ndarray         # (nlon,)

    def __post_init__(self) -> None:
        if self.df_standard.ndim != 3 or self.df_standard.shape[0] != 365:
            raise ValueError(
                f"df_standard must have shape (365, nlat, nlon); "
                f"got {self.df_standard.shape}"
            )
        if self.df_leap.ndim != 3 or self.df_leap.shape[0] != 366:
            raise ValueError(
                f"df_leap must have shape (366, nlat, nlon); "
                f"got {self.df_leap.shape}"
            )
        if self.df_standard.shape[1:] != self.df_leap.shape[1:]:
            raise ValueError(
                f"df_standard and df_leap must have the same (lat, lon) "
                f"shape; got {self.df_standard.shape[1:]} vs "
                f"{self.df_leap.shape[1:]}"
            )
        nlat = self.df_standard.shape[1]
        nlon = self.df_standard.shape[2]
        if self.source_lats_deg.shape != (nlat,):
            raise ValueError(
                f"source_lats_deg shape {self.source_lats_deg.shape} does "
                f"not match the source lat axis length {nlat}"
            )
        if self.source_lons_deg.shape != (nlon,):
            raise ValueError(
                f"source_lons_deg shape {self.source_lons_deg.shape} does "
                f"not match the source lon axis length {nlon}"
            )

    def slice_for_date(self, dt: datetime) -> np.ndarray:
        doy_0indexed = dt.timetuple().tm_yday - 1
        if calendar.isleap(dt.year):
            return self.df_leap[doy_0indexed]
        else:
            return self.df_standard[doy_0indexed]


# -----------------------------------------------------------------------------
# Summary returned by assemble_into_writer — surfaces key statistics for
# transparency logging and regression tests.
# -----------------------------------------------------------------------------

@dataclass
class AssemblyReport:
    nt: int
    nlat: int
    nlon: int
    elapsed_seconds: float
    # Per-day DF extrema and cells-outside-grid count, across all days:
    df_min: float
    df_max: float
    df_mean: float
    total_filled_cells: int
    # The Horizons-date → source-index mapping verified at day 0,
    # preserved here so a reader can confirm.
    first_date: str
    first_source_index: int
    first_source_is_leap: bool


# -----------------------------------------------------------------------------
# Pure-numpy orchestrator — the testable core.
# -----------------------------------------------------------------------------

def assemble_into_writer(
    source_provider: SourceProviderProtocol,
    ephemeris: HorizonsEphemeris,
    writer: WriterProtocol,
    subsolar_lon_deg: float = 0.0,
    fill_value: float = 1.0,
    progress_every: int = 365,
    show_progress: bool = False,
    df_plausibility_tolerance: float = 1e-6,
) -> AssemblyReport:
    """
    Rotate every day's sun-fixed slab and stream it into ``writer``.

    Parameters
    ----------
    source_provider :
        Implements :class:`SourceProviderProtocol`.  Typically a
        :class:`SunFixedSourceProvider`.
    ephemeris :
        Parsed Horizons ephemeris (see ``horizons_parser.parse``).
    writer :
        Implements :class:`WriterProtocol`.  Typically
        ``write_cesm_nc.CesmNcWriter``, but can be a test mock.
    subsolar_lon_deg :
        Subsolar longitude at the moment of each snapshot.  Defaults to 0°
        (noon UT, with substellar over the prime meridian — the convention
        established by the MATLAB upstream and recorded in the top-level
        README).
    fill_value :
        Passed through to :func:`rotate.rotate_shade_field`.  1.0 = no shade.
    progress_every :
        Log a progress line every this-many days (default 365 ≈ annual).
        Ignored when ``show_progress=True`` (use the tqdm bar instead).
    show_progress :
        If True, display a tqdm progress bar over the daily loop.  Falls
        back to ``progress_every``-style logging if ``tqdm`` is not
        installed (with a one-time warning).  Default False so library
        callers (notebooks, tests) don't get a surprise bar.
    df_plausibility_tolerance :
        A per-day slab may not exceed 1 + tolerance or go below 0 - tolerance.
        DF is an energy ratio and must lie in [0, 1].  A violation raises
        ``ValueError`` with the offending day's date, so corrupt input is
        caught early rather than silently overwriting 5 GB of disk.

    Returns
    -------
    AssemblyReport — summary statistics for the run.
    """
    nt = ephemeris.nt
    nlat = source_provider.source_lats_deg.size
    nlon = source_provider.source_lons_deg.size

    if nt == 0:
        raise ValueError("ephemeris is empty — nothing to assemble")

    # ---- Day-zero transparency log ----------------------------------------
    # Spell out exactly what gets used on day 0 so the reader can verify.
    dt0 = ephemeris.datetimes[0]
    is_leap0 = calendar.isleap(dt0.year)
    doy0 = dt0.timetuple().tm_yday - 1
    logger.info(
        "Assembly begins: %d daily snapshots (%s → %s), "
        "grid %d lat × %d lon, subsolar_lon=%.2f°",
        nt, ephemeris.date_strings[0], ephemeris.date_strings[-1],
        nlat, nlon, subsolar_lon_deg,
    )
    logger.info(
        "Day-0 source-slab mapping: Horizons date=%s → source=%s, "
        "day-of-year-0indexed=%d, δ=%+.4f°",
        ephemeris.date_strings[0],
        "df_leap" if is_leap0 else "df_standard",
        doy0,
        float(ephemeris.subsolar_lat_deg[0]),
    )

    # ---- Running statistics for the summary --------------------------------
    t_start = _time.perf_counter()
    df_min = float("+inf")
    df_max = float("-inf")
    df_sum = 0.0
    total_filled = 0

    # ---- Choose the iterator: tqdm bar or plain range --------------------
    # If show_progress=True and tqdm is importable, wrap the index range
    # with a bar.  If tqdm is missing, log a warning once and fall back to
    # the existing periodic log line.
    use_tqdm = False
    if show_progress:
        try:
            from tqdm import tqdm   # lazy import — optional dependency
            iterator = tqdm(
                range(nt),
                desc="Rotating daily shade slabs",
                unit="day",
                total=nt,
            )
            use_tqdm = True
        except ImportError:
            logger.warning(
                "show_progress=True but tqdm is not installed; falling "
                "back to log-line progress every %d days. "
                "Install with: pip install tqdm",
                progress_every,
            )
            iterator = range(nt)
    else:
        iterator = range(nt)

    # ---- Main loop --------------------------------------------------------
    for t in iterator:
        dt = ephemeris.datetimes[t]
        delta_deg = float(ephemeris.subsolar_lat_deg[t])

        slab_sun = source_provider.slice_for_date(dt)
        if slab_sun.shape != (nlat, nlon):
            raise ValueError(
                f"Source provider returned shape {slab_sun.shape} for "
                f"{dt.date()}; expected ({nlat}, {nlon})"
            )

        slab_geo, diag = rotate_shade_field(
            slab_sun,
            subsolar_lat_deg=delta_deg,
            subsolar_lon_deg=subsolar_lon_deg,
            source_lats_deg=source_provider.source_lats_deg,
            source_lons_deg=source_provider.source_lons_deg,
            fill_value=fill_value,
            return_diagnostics=True,
        )

        # Plausibility: DF is an energy ratio; it belongs to [0, 1].
        if (
            diag.output_min < -df_plausibility_tolerance
            or diag.output_max > 1.0 + df_plausibility_tolerance
        ):
            raise ValueError(
                f"Implausible DF on {dt.date()} (t={t}): "
                f"min={diag.output_min:.6f}, max={diag.output_max:.6f}. "
                f"Expected values in [0, 1]."
            )

        writer.write_day(t, slab_geo)

        # Aggregate statistics (use the diagnostics the rotate call
        # already computed — no extra passes over the array).
        if diag.output_min < df_min:
            df_min = diag.output_min
        if diag.output_max > df_max:
            df_max = diag.output_max
        df_sum += float(slab_geo.mean())
        total_filled += diag.n_filled

        # Periodic log-based progress — suppressed when tqdm is showing,
        # because the bar already tells the story.
        if (
            not use_tqdm
            and progress_every
            and (t + 1) % progress_every == 0
        ):
            elapsed = _time.perf_counter() - t_start
            rate = (t + 1) / elapsed if elapsed > 0 else 0.0
            remaining = (nt - t - 1) / rate if rate > 0 else 0.0
            logger.info(
                "  progress: %d/%d days (%.1f%%), %.1f days/s, "
                "ETA %.0fs",
                t + 1, nt, 100.0 * (t + 1) / nt, rate, remaining,
            )

    elapsed = _time.perf_counter() - t_start
    df_mean = df_sum / nt

    logger.info(
        "Assembly complete in %.2fs: DF range [%.6f, %.6f], "
        "cross-day mean of per-day means = %.6f, "
        "total filled cells = %d across all days.",
        elapsed, df_min, df_max, df_mean, total_filled,
    )

    return AssemblyReport(
        nt=nt, nlat=nlat, nlon=nlon,
        elapsed_seconds=elapsed,
        df_min=df_min, df_max=df_max, df_mean=df_mean,
        total_filled_cells=total_filled,
        first_date=ephemeris.date_strings[0],
        first_source_index=doy0,
        first_source_is_leap=is_leap0,
    )


# -----------------------------------------------------------------------------
# I/O glue — reads Matt's source file; thin wrapper around netCDF4.
# -----------------------------------------------------------------------------

# Expected variable names in Matt's NetCDF (short names, per Francis's
# solar_shade.F90 and the ArchB_df_dfl_f09_structure.nc template).
_SOURCE_VAR_DF_STANDARD = "df"
_SOURCE_VAR_DF_LEAP = "dfl"
_SOURCE_VAR_LAT = "lat"
_SOURCE_VAR_LON = "lon"


def _validate_lat_array(lat: np.ndarray, nlat_expected: int) -> bool:
    """Return True if ``lat`` looks like a real latitude axis."""
    if lat.shape != (nlat_expected,):
        return False
    if not np.all(np.isfinite(lat)):
        return False
    if lat.min() < -90.01 or lat.max() > 90.01:
        return False
    # Must be monotonic in either direction.
    diffs = np.diff(lat)
    return bool(np.all(diffs > 0) or np.all(diffs < 0))


def _validate_lon_array(lon: np.ndarray, nlon_expected: int) -> bool:
    """Return True if ``lon`` looks like a real longitude axis (0..360).

    The template file ``ArchB_df_dfl_f09_structure.nc`` writes the integer
    indices ``0, 1, 2, ..., nlon-1`` into the lon coordinate variable —
    superficially plausible (monotonic, in [0, 360]) but actually missing
    most of the globe.  The clinching diagnostic is:

        (mean Δlon) × nlon  ≈ 360°

    for any real lon-axis covering the full globe (whether or not it
    includes the wrap point).  For the f09 grid: 1.25° × 288 = 360.
    For the template: 1.0° × 288 = 288 — caught here.
    """
    if lon.shape != (nlon_expected,):
        return False
    if not np.all(np.isfinite(lon)):
        return False
    if lon.min() < -0.01 or lon.max() > 360.01:
        return False
    if not np.all(np.diff(lon) > 0):
        return False
    # The implied step times the number of points should cover the globe.
    mean_step = float(np.mean(np.diff(lon)))
    implied_coverage = mean_step * nlon_expected
    return abs(implied_coverage - 360.0) < 36.0   # within 10% of full globe


def _f09_fallback_grid(nlat: int, nlon: int) -> Tuple[np.ndarray, np.ndarray]:
    """CESM f09 grid: linspace(-90, 90, nlat) × linspace(0, 360, nlon, endpoint=False).

    This is the standard CAM FV09 grid (1° × 1.25° at 192 × 288).
    """
    return (
        np.linspace(-90.0, 90.0, nlat),
        np.linspace(0.0, 360.0, nlon, endpoint=False),
    )


def read_source_netcdf(
    path: Union[str, Path],
) -> SunFixedSourceProvider:
    """
    Read Matt's 2-year sun-fixed NetCDF and return a ``SunFixedSourceProvider``.

    Expects variables ``df`` (365, nlat, nlon), ``dfl`` (366, nlat, nlon),
    and coordinate arrays ``lat``, ``lon``.  If the lat/lon coordinate
    variables are missing or implausible, falls back to the CESM f09 grid
    with a loud warning — this is the case for the ``*_structure.nc``
    template, which stores placeholder integers in ``lon``.

    Requires ``netCDF4`` at runtime.  We import it lazily so that the rest
    of this module (the pure-numpy orchestrator) remains importable and
    testable in environments without netCDF4 installed.
    """
    try:
        import netCDF4  # type: ignore
    except ImportError as e:
        raise ImportError(
            "Reading Matt's source NetCDF requires the `netCDF4` package. "
            "Install it with `pip install netCDF4`."
        ) from e

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Source NetCDF not found: {path}")

    with netCDF4.Dataset(str(path), "r") as ds:
        missing = [
            v for v in (
                _SOURCE_VAR_DF_STANDARD, _SOURCE_VAR_DF_LEAP,
            ) if v not in ds.variables
        ]
        if missing:
            raise KeyError(
                f"Source NetCDF {path} is missing expected variables: "
                f"{missing}. Found variables: {list(ds.variables.keys())}"
            )

        df_standard = np.asarray(ds[_SOURCE_VAR_DF_STANDARD][:], dtype=np.float64)
        df_leap = np.asarray(ds[_SOURCE_VAR_DF_LEAP][:], dtype=np.float64)

        # Handle (lat, lon, time) axis ordering that Matt's MATLAB ncwrite
        # may produce.  We want (time, lat, lon) for Python use.  If the
        # leading axis is 365/366, we're good; otherwise we need to move it.
        df_standard = _maybe_reorder_to_time_first(df_standard, expect_time=365)
        df_leap = _maybe_reorder_to_time_first(df_leap, expect_time=366)

        nlat = df_standard.shape[1]
        nlon = df_standard.shape[2]

        lat = None
        lon = None
        if _SOURCE_VAR_LAT in ds.variables:
            lat_candidate = np.asarray(ds[_SOURCE_VAR_LAT][:], dtype=np.float64)
            if _validate_lat_array(lat_candidate, nlat):
                lat = lat_candidate
        if _SOURCE_VAR_LON in ds.variables:
            lon_candidate = np.asarray(ds[_SOURCE_VAR_LON][:], dtype=np.float64)
            if _validate_lon_array(lon_candidate, nlon):
                lon = lon_candidate

    if lat is None or lon is None:
        f09_lat, f09_lon = _f09_fallback_grid(nlat, nlon)
        logger.warning(
            "Source NetCDF %s has missing or implausible lat/lon coordinate "
            "variables; falling back to CESM f09 grid "
            "(lat: linspace(-90, 90, %d); lon: linspace(0, 360, %d, endpoint=False)).",
            path, nlat, nlon,
        )
        if lat is None:
            lat = f09_lat
        if lon is None:
            lon = f09_lon

    # Ensure lat is south-to-north (monotonically increasing); flip data if not.
    if lat[0] > lat[-1]:
        logger.info(
            "Source lat axis is decreasing (north-to-south); flipping "
            "both lat and the lat axis of df/dfl to be south-to-north.",
        )
        lat = lat[::-1].copy()
        df_standard = df_standard[:, ::-1, :].copy()
        df_leap = df_leap[:, ::-1, :].copy()

    logger.info(
        "Loaded source NetCDF %s: df %s, dfl %s, lat[%d..%d]=[%.3f, %.3f], "
        "lon[%d..%d]=[%.3f, %.3f]",
        path, df_standard.shape, df_leap.shape,
        0, nlat - 1, lat[0], lat[-1],
        0, nlon - 1, lon[0], lon[-1],
    )

    return SunFixedSourceProvider(
        df_standard=df_standard, df_leap=df_leap,
        source_lats_deg=lat, source_lons_deg=lon,
    )


def _maybe_reorder_to_time_first(
    arr: np.ndarray, expect_time: int,
) -> np.ndarray:
    """
    Return ``arr`` with the time axis first.  Accepts any 3-D array with one
    axis of length ``expect_time`` (365 or 366) — we move it to position 0
    and leave the remaining two axes in their original order.

    Matt's MATLAB ``ncwrite`` stores in Fortran order so a variable declared
    ``(lon=288, lat=192, time=365)`` reads in Python as
    ``(time=365, lat=192, lon=288)`` — already time-first — but we
    nevertheless check, so this will handle both conventions correctly.
    """
    if arr.ndim != 3:
        raise ValueError(f"Expected 3-D source array; got shape {arr.shape}")
    if arr.shape[0] == expect_time:
        return arr
    # Find the time axis.
    axes = [i for i, s in enumerate(arr.shape) if s == expect_time]
    if len(axes) != 1:
        raise ValueError(
            f"Cannot unambiguously identify the time axis in shape "
            f"{arr.shape}: no axis (or multiple axes) have length "
            f"{expect_time}."
        )
    t_axis = axes[0]
    logger.warning(
        "Source array shape %s has time axis at position %d; moving it to "
        "position 0 for consistent (time, lat, lon) ordering.",
        arr.shape, t_axis,
    )
    return np.moveaxis(arr, t_axis, 0)


# -----------------------------------------------------------------------------
# Top-level convenience: from input paths to output file.
# -----------------------------------------------------------------------------

def assemble_from_files(
    source_nc_path: Union[str, Path],
    horizons_txt_path: Union[str, Path],
    output_nc_path: Union[str, Path],
    subsolar_lon_deg: float = 0.0,
    fill_value: float = 1.0,
    show_progress: bool = False,
) -> AssemblyReport:
    """
    Read, assemble, write — the one-call "run it" convenience.

    Uses :func:`read_source_netcdf` (requires netCDF4) for input and
    ``write_cesm_nc.CesmNcWriter`` for output.

    The two input paths (source NetCDF, Horizons ephemeris) are forwarded
    to the writer so they land in the output file's ``history`` global
    attribute — anyone opening the resulting ``.nc`` with ``ncdump -h``
    can trace it back to the exact inputs that produced it.
    """
    source = read_source_netcdf(source_nc_path)
    ephemeris = horizons_parser.parse(horizons_txt_path)

    # Lazy import so the rest of assemble.py can be used without the
    # writer being present (useful for tests and partial environments).
    from .write_cesm_nc import CesmNcWriter

    with CesmNcWriter(
        path=output_nc_path,
        source_lats_deg=source.source_lats_deg,
        source_lons_deg=source.source_lons_deg,
        ephemeris=ephemeris,
        source_nc_path=source_nc_path,
        horizons_path=horizons_txt_path,
    ) as writer:
        return assemble_into_writer(
            source_provider=source,
            ephemeris=ephemeris,
            writer=writer,
            subsolar_lon_deg=subsolar_lon_deg,
            fill_value=fill_value,
            show_progress=show_progress,
        )


# -----------------------------------------------------------------------------
# Script entry point — one-off runs.
# -----------------------------------------------------------------------------

def _main(argv: Optional[List[str]] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(
        description=(
            "Rotate Matt's sun-fixed DF NetCDF into geographic coordinates, "
            "one day per Horizons ephemeris entry."
        )
    )
    ap.add_argument("source_nc", help="Path to Matt's 2-year DF NetCDF.")
    ap.add_argument("horizons_txt", help="Path to the JPL Horizons ephemeris.")
    ap.add_argument("output_nc", help="Path to write the 30-year DF NetCDF.")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    report = assemble_from_files(
        source_nc_path=args.source_nc,
        horizons_txt_path=args.horizons_txt,
        output_nc_path=args.output_nc,
        show_progress=True,   # one-off script: show the tqdm bar
    )
    logger.info("Report: %s", report)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
