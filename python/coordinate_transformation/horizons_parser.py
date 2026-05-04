"""
horizons_parser
===============

Parse a JPL Horizons ephemeris text file and return the Sun's declination
(subsolar latitude) for each day.

Background
----------
The coordinate transformation from the PSF sun-fixed frame into geographic
(latitude, longitude) coordinates requires, for every simulated day, the Sun's
declination δ at the moment the model computes the shade pattern (00:00 UT by
convention). JPL's Horizons system (https://ssd.jpl.nasa.gov/horizons/) is the
standard source of solar ephemeris data. The file parsed here was generated
with:

    Target body      : Sun (10)
    Center body      : Earth (399)
    Center-site      : GEOCENTRIC
    Start/Stop       : 2035-01-01 through 2065-12-31
    Step-size        : 1440 minutes (i.e. one step per day at 00:00 UT)
    Output           : RA and Declination, sexagesimal (HH MM SS, ±DD MM SS)

The ``declination`` column is the equatorial declination of the Sun in the
geocentric frame — equivalently, the *subsolar latitude*: the geographic
latitude at which the Sun is directly overhead at that moment. That is exactly
the quantity Illeana's rotation routine needs.

File layout
-----------
The data section is bracketed by the Horizons markers ``$$SOE`` (Start Of
Ephemeris) and ``$$EOE`` (End Of Ephemeris). Each line between them has the
form::

     2035-Jan-01 00:00     18 43 02.53 -23 03 57.6
    └──┬──────┘ └──┬┘   └────┬──────┘ └────┬──────┘
       date       time        RA (HMS)     Dec (±DMS)

Only the date and the declination are needed; RA is ignored.

Sign convention for declination
-------------------------------
Horizons stores declination with the sign on the degree field. A value of
``-23 03 57.6`` means::

    δ = -(23° + 3'/60 + 57.6"/3600) = -23.0660°

**not**::

    δ = -23° + 3' + 57.6"   (which would be -22.9340°)

Arc-minute and arc-second fields are always non-negative magnitude components.
The parser below handles this correctly, including the edge case where the
degree field reads ``-00`` (a small negative declination between 0° and -1°).

Output
------
:func:`parse` returns a :class:`HorizonsEphemeris` with:

* ``nt``               — number of days parsed
* ``date_strings``     — ISO ``YYYY-MM-DD`` strings, length ``nt``
* ``datetimes``        — Python ``datetime`` objects at 00:00 UT, length ``nt``
* ``subsolar_lat_deg`` — signed decimal degrees, ``np.float64`` array of
  length ``nt``. Positive = northern hemisphere subsolar point.

On parse, the first entry is logged at INFO level — a transparency requirement
per Illeana's spec ("make it write an example of what it is taking").
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Union

import numpy as np


logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

# Horizons uses three-letter English month abbreviations. We map them
# explicitly here rather than relying on datetime.strptime('%b', ...), which
# is locale-dependent and would silently misparse on a non-English system.
_MONTH_ABBR_TO_NUM = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
}

_SOE_MARKER = "$$SOE"
_EOE_MARKER = "$$EOE"


# -----------------------------------------------------------------------------
# Return type
# -----------------------------------------------------------------------------

@dataclass
class HorizonsEphemeris:
    """Parsed contents of a Horizons solar-ephemeris text file."""

    nt: int
    date_strings: List[str]
    datetimes: List[datetime]
    subsolar_lat_deg: np.ndarray  # shape (nt,), float64, signed

    def __repr__(self) -> str:  # short, for interactive use
        if self.nt == 0:
            return "HorizonsEphemeris(empty)"
        return (
            f"HorizonsEphemeris(nt={self.nt}, "
            f"first={self.date_strings[0]} δ={self.subsolar_lat_deg[0]:+.4f}°, "
            f"last={self.date_strings[-1]} δ={self.subsolar_lat_deg[-1]:+.4f}°)"
        )


# -----------------------------------------------------------------------------
# Sub-parsers (exposed so tests can exercise them in isolation)
# -----------------------------------------------------------------------------

def parse_declination_dms(deg_field: str, arcmin_field: str, arcsec_field: str) -> float:
    """
    Convert a Horizons ±DD MM SS.s declination triplet into signed decimal
    degrees.

    Parameters
    ----------
    deg_field   : e.g. ``"-23"``, ``"+00"``, ``"-00"``, ``"23"``
    arcmin_field: e.g. ``"03"``; always a non-negative magnitude
    arcsec_field: e.g. ``"57.6"``; always a non-negative magnitude

    Returns
    -------
    δ in decimal degrees (float), signed.

    Notes
    -----
    The sign on the degree field applies to the whole sexagesimal value — see
    the module docstring. ``-00 30 00.0`` therefore parses to ``-0.5000°``,
    not ``+0.5000°``.
    """
    # Detect explicit negative sign. Any '+' or no-sign is treated as positive.
    # abs(int(...)) then handles both '+00', '-00', '23', and '-23' uniformly.
    is_negative = deg_field.lstrip().startswith("-")
    deg_magnitude = abs(int(deg_field))
    arcmin = int(arcmin_field)
    arcsec = float(arcsec_field)

    if arcmin < 0 or arcsec < 0:
        raise ValueError(
            f"Arcmin and arcsec fields must be non-negative in Horizons format; "
            f"got arcmin={arcmin_field!r}, arcsec={arcsec_field!r}"
        )

    magnitude = deg_magnitude + arcmin / 60.0 + arcsec / 3600.0
    return -magnitude if is_negative else magnitude


def parse_date_token(token: str) -> datetime:
    """
    Parse a Horizons date token (``YYYY-Mon-DD``, e.g. ``"2035-Jan-01"``) into
    a Python ``datetime`` at 00:00 UT.
    """
    parts = token.split("-")
    if len(parts) != 3:
        raise ValueError(f"Expected YYYY-Mon-DD, got {token!r}")
    year_str, mon_str, day_str = parts
    if mon_str not in _MONTH_ABBR_TO_NUM:
        raise ValueError(
            f"Unknown month abbreviation {mon_str!r} in date token {token!r}"
        )
    return datetime(int(year_str), _MONTH_ABBR_TO_NUM[mon_str], int(day_str))


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

def parse(
    path: Union[str, Path],
    log_first_entry: bool = True,
) -> HorizonsEphemeris:
    """
    Read a Horizons ephemeris file and return the parsed ephemeris.

    Parameters
    ----------
    path : path to the Horizons ``.txt`` file.
    log_first_entry : if True (default), log the first parsed row at INFO
        level so the user can confirm the parser sees what they expect.
        Illeana's spec asks for this: "make it write an example of what it is
        taking."

    Returns
    -------
    HorizonsEphemeris
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Horizons file not found: {path}")

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    # Locate the SOE/EOE markers. Horizons writes them on their own lines.
    soe_idx, eoe_idx = None, None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped == _SOE_MARKER and soe_idx is None:
            soe_idx = i
        elif stripped == _EOE_MARKER:
            eoe_idx = i
            break

    if soe_idx is None:
        raise ValueError(f"{_SOE_MARKER} marker not found in {path}")
    if eoe_idx is None:
        raise ValueError(f"{_EOE_MARKER} marker not found in {path}")
    if eoe_idx <= soe_idx + 1:
        raise ValueError(
            f"No ephemeris rows between {_SOE_MARKER} (line {soe_idx + 1}) "
            f"and {_EOE_MARKER} (line {eoe_idx + 1}) in {path}"
        )

    date_strings: List[str] = []
    datetimes: List[datetime] = []
    decs: List[float] = []

    # Parse every line strictly between SOE and EOE.
    for idx in range(soe_idx + 1, eoe_idx):
        raw = lines[idx]
        tokens = raw.split()

        # A well-formed ephemeris row has at least:
        #   date, time, RA_h, RA_m, RA_s, Dec_deg, Dec_arcmin, Dec_arcsec
        # → 8 tokens. We parse by position, not regex, because Horizons
        # delimits with runs of whitespace (multiple spaces) and tokens are
        # unambiguous once split.
        if len(tokens) < 8:
            raise ValueError(
                f"Unexpected ephemeris row format at file line {idx + 1} in "
                f"{path!s}:\n  {raw.rstrip()!r}\n"
                f"Expected at least 8 whitespace-separated tokens "
                f"(date, time, RA HMS, Dec DMS); got {len(tokens)}."
            )

        date_token = tokens[0]
        dt = parse_date_token(date_token)

        # tokens[1] is the time-of-day field (e.g. '00:00') — redundant here
        # because we requested a daily step at 00:00 UT. We skip it.
        # tokens[2..4] are RA h/m/s — ignored for subsolar latitude.

        dec_deg = parse_declination_dms(tokens[5], tokens[6], tokens[7])

        date_strings.append(dt.strftime("%Y-%m-%d"))
        datetimes.append(dt)
        decs.append(dec_deg)

    nt = len(date_strings)
    result = HorizonsEphemeris(
        nt=nt,
        date_strings=date_strings,
        datetimes=datetimes,
        subsolar_lat_deg=np.asarray(decs, dtype=np.float64),
    )

    if log_first_entry and nt > 0:
        logger.info(
            "Parsed %d daily ephemeris entries from %s. "
            "Example of what it is taking: "
            "date=%s (datetime=%s), subsolar_lat=%+.4f°",
            nt,
            path,
            result.date_strings[0],
            result.datetimes[0].isoformat(),
            float(result.subsolar_lat_deg[0]),
        )

    return result


# -----------------------------------------------------------------------------
# Command-line entry — for quick standalone inspection.
# -----------------------------------------------------------------------------

def _main() -> int:
    import argparse

    ap = argparse.ArgumentParser(
        description="Parse a JPL Horizons solar-ephemeris text file and "
                    "print summary statistics of the subsolar latitude."
    )
    ap.add_argument("path", help="Path to the Horizons .txt file.")
    ap.add_argument(
        "--show", type=int, default=3,
        help="Number of first/last rows to display (default: 3).",
    )
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    eph = parse(args.path)

    print(f"Parsed {eph.nt} daily entries.")
    print("First rows:")
    for i in range(min(args.show, eph.nt)):
        print(f"  {eph.date_strings[i]}   δ = {eph.subsolar_lat_deg[i]:+.4f}°")
    print("Last rows:")
    for i in range(max(0, eph.nt - args.show), eph.nt):
        print(f"  {eph.date_strings[i]}   δ = {eph.subsolar_lat_deg[i]:+.4f}°")

    lat = eph.subsolar_lat_deg
    print(
        f"Summary: δ range [{lat.min():+.4f}°, {lat.max():+.4f}°]; "
        f"mean {lat.mean():+.4f}°; "
        f"±23.44° expected (obliquity of the ecliptic)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
