"""
write_cesm_nc
=============

Write the geographic-frame dimming-factor field into a CF-compliant NetCDF4
file.  The writer is used by :mod:`assemble` via the ``WriterProtocol``
interface defined there.

File schema
-----------
::

    dimensions:
        time       = <nt>                  # unlimited-like, but fixed here
        latitude   = <nlat>                # e.g. 192 for f09
        longitude  = <nlon>                # e.g. 288 for f09
        string_len = 10                    # len("YYYY-MM-DD")

    variables:
        double time(time)
            units         = "days since YYYY-MM-DD 12:00:00"
            calendar      = "gregorian"
            standard_name = "time"
            long_name     = "time (12:00 UT on each modelled day)"
            axis          = "T"

        double latitude(latitude)
            units         = "degrees_north"
            standard_name = "latitude"
            long_name     = "latitude"
            axis          = "Y"

        double longitude(longitude)
            units         = "degrees_east"
            standard_name = "longitude"
            long_name     = "longitude"
            axis          = "X"

        char date(time, string_len)
            long_name     = "calendar date (ISO 8601)"

        double DF(time, latitude, longitude)
            long_name     = "dimming factor"
            units         = "1"                    # dimensionless
            valid_range   = (0.0, 1.0)
            _FillValue    = 9.96921e36             # NetCDF default f64 fill
            comment       = "Fraction of surface solar irradiance passing "
                            "the sunshade constellation. 1.0 = no shade; "
                            "0.0 = total occultation."

    global attributes:
        Conventions = "CF-1.8"
        title, institution, source, references, history, comment

Transparency
------------
The point of this file is that **every number in it should be traceable**.
The global ``history`` attribute carries the inputs and UTC creation time;
``source`` names the rotation algorithm module path; ``references`` points
at the project repository.  A scientist opening the output with ``ncdump
-h`` should see enough to reconstruct what produced it.

Atomic writes
-------------
To avoid leaving a half-written 5 GB file on disk after a crash, the
writer writes to ``<output>.partial`` and renames to the final path only
when :meth:`close` succeeds.  If the writer is used as a context manager
and an exception is raised inside the ``with`` block, the partial file is
removed.  The final file therefore either exists fully or not at all.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import numpy as np

from .horizons_parser import HorizonsEphemeris


logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Constants — the schema is pinned in one place for easy review.
# -----------------------------------------------------------------------------

DIM_TIME = "time"
DIM_LAT = "latitude"
DIM_LON = "longitude"
DIM_STRING_LEN = "string_len"
STRING_LEN = 10                        # len("YYYY-MM-DD")

VAR_TIME = "time"
VAR_LAT = "latitude"
VAR_LON = "longitude"
VAR_DATE = "date"
VAR_DF = "DF"

DF_FILL_VALUE = 9.96921e36             # NetCDF's default f64 _FillValue
DF_COMPRESSION_LEVEL = 4               # 0=off, 1=fastest, 9=smallest

CF_CONVENTIONS = "CF-1.8"
REPO_URL = (
    "https://github.com/PlanetarySunshadeInstitute/"
    "sunshade-irradiance-model"
)


# -----------------------------------------------------------------------------
# Pure helpers — no netCDF4 dependency
# -----------------------------------------------------------------------------

def encode_cf_time(
    datetimes: Sequence[datetime],
    reference: Optional[datetime] = None,
) -> Tuple[np.ndarray, str]:
    """
    Convert a sequence of ``datetime`` objects into CF "days since ..."
    numeric values.

    Parameters
    ----------
    datetimes : sequence of datetime
        The timestamps to encode.  Naïve (timezone-unaware) is expected,
        consistent with how the Horizons parser returns them.
    reference : datetime, optional
        The zero point.  Defaults to ``datetimes[0]`` truncated to the
        start of that day at 00:00:00.

    Returns
    -------
    values : ndarray of float64, shape (len(datetimes),)
        Days since the reference.  Float so sub-daily resolution is
        representable if ever needed.  In production the writer call
        site (:meth:`CesmNcWriter._write_header`) shifts the Horizons
        daily timestamps to noon UT and passes a noon-UT reference, so
        the written values are integer-valued days.
    units_str : str
        The CF ``units`` attribute string, e.g. "days since 2035-01-01 00:00:00".
    """
    if len(datetimes) == 0:
        raise ValueError("encode_cf_time requires at least one timestamp")
    if reference is None:
        d0 = datetimes[0]
        reference = datetime(d0.year, d0.month, d0.day)
    one_day = 86400.0
    # Use total_seconds on the timedelta so sub-daily resolution works.
    values = np.asarray(
        [(dt - reference).total_seconds() / one_day for dt in datetimes],
        dtype=np.float64,
    )
    units_str = f"days since {reference.strftime('%Y-%m-%d %H:%M:%S')}"
    return values, units_str


def build_global_attrs(
    ephemeris: HorizonsEphemeris,
    source_nc_path: Optional[Union[str, Path]] = None,
    horizons_path: Optional[Union[str, Path]] = None,
    extra_history: str = "",
) -> dict:
    """
    Build the global-attribute dictionary written at the top of the file.

    Every field here is fixed-at-time-of-write so the file is completely
    self-describing (``ncdump -h`` gives a reader everything they need to
    understand what produced the data).
    """
    created_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    history_lines = [
        f"{created_utc}: created by coordinate_transformation.assemble"
    ]
    if source_nc_path is not None:
        history_lines.append(f"  source NetCDF: {source_nc_path}")
    if horizons_path is not None:
        history_lines.append(f"  ephemeris: {horizons_path}")
    history_lines.append(
        f"  period: {ephemeris.date_strings[0]} through "
        f"{ephemeris.date_strings[-1]} ({ephemeris.nt} daily snapshots)"
    )
    if extra_history:
        history_lines.append(extra_history)

    return {
        "Conventions": CF_CONVENTIONS,
        "title": (
            "Planetary Sunshade dimming factor in geographic coordinates"
        ),
        "institution": "Planetary Sunshade Institute",
        "source": (
            "PSF sun-fixed irradiance model rotated into geographic "
            "coordinates using Sun's declination from JPL Horizons; "
            "see module coordinate_transformation.rotate."
        ),
        "references": REPO_URL,
        "history": "\n".join(history_lines),
        "comment": (
            "DF is the fraction of surface solar irradiance that passes "
            "the sunshade constellation, computed once per day at noon UT "
            "with the Sun's substellar point at geographic longitude 0°. "
            "Intended to be multiplied, externally to CESM, against the "
            "CMIP6 SolarForcing TSI/SSI file to produce a modulated "
            "spectral irradiance input for climate simulation."
        ),
    }


def resolve_atomic_paths(
    final_path: Union[str, Path],
) -> Tuple[Path, Path]:
    """
    Return ``(final_path, partial_path)``.  The partial path is the final
    path with a ``.partial`` suffix appended; writers stream into this and
    rename on successful close.
    """
    final_path = Path(final_path)
    partial_path = final_path.with_name(final_path.name + ".partial")
    return final_path, partial_path


# -----------------------------------------------------------------------------
# Writer
# -----------------------------------------------------------------------------

@dataclass
class WriterOptions:
    """Tunable knobs for the writer, bundled so the constructor stays
    readable and so tests can override defaults cleanly."""

    compress: bool = True
    compress_level: int = DF_COMPRESSION_LEVEL
    nc_format: str = "NETCDF4"
    remove_partial_on_error: bool = True


class CesmNcWriter:
    """
    Streams the rotated dimming-factor field into a CF-compliant NetCDF.

    Usage
    -----
    ::

        with CesmNcWriter(
            path="df_2035_2065.nc",
            source_lats_deg=lats, source_lons_deg=lons,
            ephemeris=ephem,
        ) as writer:
            for t in range(nt):
                writer.write_day(t, rotated_slab_2d)

    The constructor creates the target file immediately (at
    ``<path>.partial``), writes the full header (dimensions, coordinates,
    date strings, attributes), and pre-allocates ``DF`` as a chunked
    (1, nlat, nlon) variable — so ``write_day`` is a single hyperslab
    write per day.
    """

    def __init__(
        self,
        path: Union[str, Path],
        source_lats_deg: np.ndarray,
        source_lons_deg: np.ndarray,
        ephemeris: HorizonsEphemeris,
        source_nc_path: Optional[Union[str, Path]] = None,
        horizons_path: Optional[Union[str, Path]] = None,
        options: Optional[WriterOptions] = None,
    ) -> None:
        try:
            import netCDF4  # type: ignore
        except ImportError as e:
            raise ImportError(
                "CesmNcWriter requires the `netCDF4` package. "
                "Install it with `pip install netCDF4`."
            ) from e

        self._nc = netCDF4  # keep a handle on the module for close/teardown
        self._options = options or WriterOptions()
        self._final_path, self._partial_path = resolve_atomic_paths(path)
        self._ephemeris = ephemeris
        self._closed = False

        # Validate inputs before opening the file — "fail fast" policy
        # means we never leave a stale .partial from rejected inputs.
        lats = np.asarray(source_lats_deg, dtype=np.float64)
        lons = np.asarray(source_lons_deg, dtype=np.float64)
        self._validate_inputs(lats, lons, ephemeris)

        self._nlat = lats.size
        self._nlon = lons.size
        self._nt = ephemeris.nt

        # Open the .partial file and write the header.
        self._ds = netCDF4.Dataset(
            str(self._partial_path),
            mode="w",
            format=self._options.nc_format,
        )
        try:
            self._write_header(lats, lons, ephemeris,
                               source_nc_path=source_nc_path,
                               horizons_path=horizons_path)
        except Exception:
            # If header writing fails, close and remove the partial.
            self._ds.close()
            if self._options.remove_partial_on_error:
                try:
                    self._partial_path.unlink()
                except FileNotFoundError:
                    pass
            raise

        logger.info(
            "Opened %s for streaming DF writes: %d time × %d lat × %d lon, "
            "compress=%s (level=%d)",
            self._partial_path, self._nt, self._nlat, self._nlon,
            self._options.compress, self._options.compress_level,
        )

    # ---- Validation ------------------------------------------------------

    @staticmethod
    def _validate_inputs(
        lats: np.ndarray, lons: np.ndarray, ephemeris: HorizonsEphemeris,
    ) -> None:
        if lats.ndim != 1:
            raise ValueError(f"source_lats_deg must be 1-D; got shape {lats.shape}")
        if lons.ndim != 1:
            raise ValueError(f"source_lons_deg must be 1-D; got shape {lons.shape}")
        if lats.size < 2 or lons.size < 2:
            raise ValueError("source grid needs at least 2 points per axis")
        if not np.all(np.diff(lats) > 0):
            raise ValueError("source_lats_deg must be strictly increasing")
        if not np.all(np.diff(lons) > 0):
            raise ValueError("source_lons_deg must be strictly increasing")
        if lats.min() < -90.01 or lats.max() > 90.01:
            raise ValueError(
                f"source_lats_deg out of range [-90, 90]: "
                f"[{lats.min()}, {lats.max()}]"
            )
        if lons.min() < -0.01 or lons.max() > 360.01:
            raise ValueError(
                f"source_lons_deg out of range [0, 360]: "
                f"[{lons.min()}, {lons.max()}]"
            )
        if ephemeris.nt == 0:
            raise ValueError("ephemeris is empty; nothing to write")
        if len(ephemeris.datetimes) != ephemeris.nt:
            raise ValueError(
                f"ephemeris.datetimes length {len(ephemeris.datetimes)} "
                f"does not match nt={ephemeris.nt}"
            )

    # ---- Header -----------------------------------------------------------

    def _write_header(
        self,
        lats: np.ndarray,
        lons: np.ndarray,
        ephemeris: HorizonsEphemeris,
        source_nc_path: Optional[Union[str, Path]],
        horizons_path: Optional[Union[str, Path]],
    ) -> None:
        ds = self._ds

        # --- Dimensions ---
        ds.createDimension(DIM_TIME, ephemeris.nt)
        ds.createDimension(DIM_LAT, lats.size)
        ds.createDimension(DIM_LON, lons.size)
        ds.createDimension(DIM_STRING_LEN, STRING_LEN)

        # --- Coordinate variables ---
        # The Horizons ephemeris is sampled at 00:00 UT, but Matt's MATLAB
        # irradiance model evaluates DF at noon UT each day (project
        # "Option A" convention — see top-level Matlab/README.md).  To keep
        # the written file internally consistent with the global-attribute
        # comment ("computed once per day at noon UT ..."), we relabel each
        # daily slab to its noon-UT timestamp before encoding, and anchor
        # the CF reference epoch at noon of the first day.  Numeric values
        # remain integer-valued days (0, 1, 2, ...).
        d0 = ephemeris.datetimes[0]
        noon_reference = datetime(d0.year, d0.month, d0.day, 12, 0, 0)
        noon_datetimes = [
            datetime(dt.year, dt.month, dt.day, 12, 0, 0)
            for dt in ephemeris.datetimes
        ]
        time_values, time_units = encode_cf_time(
            noon_datetimes, reference=noon_reference,
        )
        v_time = ds.createVariable(VAR_TIME, "f8", (DIM_TIME,))
        v_time.units = time_units
        v_time.calendar = "gregorian"
        v_time.standard_name = "time"
        v_time.long_name = "time (12:00 UT on each modelled day)"
        v_time.axis = "T"
        v_time[:] = time_values

        v_lat = ds.createVariable(VAR_LAT, "f8", (DIM_LAT,))
        v_lat.units = "degrees_north"
        v_lat.standard_name = "latitude"
        v_lat.long_name = "latitude"
        v_lat.axis = "Y"
        v_lat[:] = lats

        v_lon = ds.createVariable(VAR_LON, "f8", (DIM_LON,))
        v_lon.units = "degrees_east"
        v_lon.standard_name = "longitude"
        v_lon.long_name = "longitude"
        v_lon.axis = "X"
        v_lon[:] = lons

        # --- Date strings as char(time, string_len) ---
        v_date = ds.createVariable(
            VAR_DATE, "S1", (DIM_TIME, DIM_STRING_LEN),
        )
        v_date.long_name = "calendar date (ISO 8601)"
        # netCDF4 helper to pack Python strings into (N, STRING_LEN) chars.
        date_chars = self._nc.stringtochar(
            np.asarray(ephemeris.date_strings, dtype=f"S{STRING_LEN}"),
        )
        v_date[:, :] = date_chars

        # --- DF(time, lat, lon) — pre-allocated, chunked, compressed ---
        df_kwargs = dict(
            dimensions=(DIM_TIME, DIM_LAT, DIM_LON),
            fill_value=DF_FILL_VALUE,
            chunksizes=(1, lats.size, lons.size),
        )
        if self._options.compress:
            df_kwargs.update(
                zlib=True, complevel=self._options.compress_level,
                shuffle=True,
            )
        v_df = ds.createVariable(VAR_DF, "f8", **df_kwargs)
        v_df.long_name = "dimming factor"
        v_df.units = "1"
        v_df.valid_range = np.array([0.0, 1.0], dtype=np.float64)
        v_df.comment = (
            "Fraction of surface solar irradiance passing the sunshade "
            "constellation. 1.0 = no shade; 0.0 = total occultation."
        )
        self._v_df = v_df

        # --- Global attributes ---
        attrs = build_global_attrs(
            ephemeris,
            source_nc_path=source_nc_path,
            horizons_path=horizons_path,
        )
        for k, v in attrs.items():
            ds.setncattr(k, v)

        # Flush the header so that ncdump against the .partial file works
        # even mid-run.
        ds.sync()

    # ---- Streaming writes ------------------------------------------------

    def write_day(self, t_idx: int, df_slab: np.ndarray) -> None:
        """
        Write the rotated 2-D slab for output time index ``t_idx``.

        Parameters
        ----------
        t_idx : int, 0-based
        df_slab : ndarray of shape (nlat, nlon), float64
        """
        if self._closed:
            raise RuntimeError("Writer is closed; no more writes accepted.")
        if not (0 <= t_idx < self._nt):
            raise IndexError(
                f"t_idx={t_idx} out of range for nt={self._nt}"
            )
        if df_slab.shape != (self._nlat, self._nlon):
            raise ValueError(
                f"df_slab shape {df_slab.shape} does not match "
                f"expected ({self._nlat}, {self._nlon})"
            )
        # Hyperslab write into the pre-allocated variable.
        self._v_df[t_idx, :, :] = df_slab

    # ---- Lifecycle -------------------------------------------------------

    def close(self) -> None:
        """Finalise the file: close the dataset, rename partial → final."""
        if self._closed:
            return
        try:
            self._ds.close()
        finally:
            self._closed = True
        # Only rename on successful close; leave .partial behind on error
        # for debugging.
        os.replace(self._partial_path, self._final_path)
        logger.info("Closed and finalised %s", self._final_path)

    def abort(self) -> None:
        """Close without finalising, removing the partial file."""
        if self._closed:
            return
        try:
            self._ds.close()
        finally:
            self._closed = True
        if self._options.remove_partial_on_error:
            try:
                self._partial_path.unlink()
            except FileNotFoundError:
                pass

    # ---- Context-manager protocol ---------------------------------------

    def __enter__(self) -> "CesmNcWriter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if exc_type is None:
            self.close()
        else:
            logger.error(
                "Exception inside CesmNcWriter context (%s: %s); aborting "
                "and removing partial file %s",
                exc_type.__name__, exc_val, self._partial_path,
            )
            self.abort()
        # Do not suppress the exception.
        return None

    def __del__(self) -> None:
        # Best-effort cleanup if someone forgets to use a with-block.
        if not self._closed:
            try:
                self.abort()
            except Exception:
                pass

    # ---- Introspection (useful in tests) --------------------------------

    @property
    def final_path(self) -> Path:
        return self._final_path

    @property
    def partial_path(self) -> Path:
        return self._partial_path

    @property
    def nt(self) -> int:
        return self._nt
