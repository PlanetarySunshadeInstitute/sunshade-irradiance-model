"""
Tests for ``horizons_parser``.

Two layers:

1. **Unit tests** on the sub-parsers :func:`parse_declination_dms` and
   :func:`parse_date_token`. These pin down the sign convention and the
   locale-independent month handling.

2. **Integration tests** against the actual Horizons file shipped in
   ``Coordinate Transformation/horizons_Sun_RA-DEC.txt``. These verify that
   the parser reads the real 2035-2065 ephemeris correctly — not just a
   synthetic fixture — and pick out astronomically meaningful rows
   (solstice, near-equinox sign crossing, endpoints) so that a scientist
   reading the tests can see exactly why each expected value is what it is.

Run standalone::

    python -m unittest python.coordinate_transformation.tests.test_horizons_parser

or::

    python python/coordinate_transformation/tests/test_horizons_parser.py

Uses only the Python standard library (``unittest``) plus NumPy, so no
pytest is required.
"""

from __future__ import annotations

import io
import logging
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import numpy as np

# Make the package importable whether the tests are run from the repo root
# (python -m unittest ...) or directly (python test_horizons_parser.py).
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent                 # .../coordinate_transformation
_PYTHON_DIR = _PKG_DIR.parent           # .../python
_REPO_ROOT = _PYTHON_DIR.parent         # .../Matlab
for p in (_REPO_ROOT, _PYTHON_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from coordinate_transformation.horizons_parser import (  # noqa: E402
    HorizonsEphemeris,
    parse,
    parse_date_token,
    parse_declination_dms,
)


# The real ephemeris file that ships with the project.
HORIZONS_TXT = (
    _REPO_ROOT / "Coordinate Transformation" / "horizons_Sun_RA-DEC.txt"
)


# =====================================================================
# Sub-parser: parse_declination_dms
# =====================================================================

class TestParseDeclinationDms(unittest.TestCase):
    """Pin down the Horizons sign convention on ±DD MM SS.s triplets."""

    # Tolerance: we're converting to decimal degrees, and the arcsec input
    # has 0.1" resolution → ~3e-5 deg. 1e-6 deg is comfortably tighter.
    TOL = 1e-6

    def test_negative_degree_applies_to_whole_value(self):
        # The canonical case that trips people up:
        # "-23 03 57.6" means -23.0660°, not -(23°) + 3' + 57.6" = -22.9340°.
        got = parse_declination_dms("-23", "03", "57.6")
        expected = -(23 + 3 / 60 + 57.6 / 3600)
        self.assertAlmostEqual(got, expected, delta=self.TOL)
        self.assertLess(got, -23.0)  # definitely more negative than -23

    def test_positive_degree_no_sign(self):
        # Horizons sometimes omits the leading '+', sometimes includes it.
        self.assertAlmostEqual(
            parse_declination_dms("23", "25", "51.7"),
            23 + 25 / 60 + 51.7 / 3600,
            delta=self.TOL,
        )

    def test_positive_degree_with_explicit_plus(self):
        self.assertAlmostEqual(
            parse_declination_dms("+23", "25", "51.7"),
            23 + 25 / 60 + 51.7 / 3600,
            delta=self.TOL,
        )

    def test_minus_zero_is_small_negative(self):
        # The edge case where δ is between 0 and -1°: degree field reads
        # "-00". The int('-00') == 0 trick alone would lose the sign.
        got = parse_declination_dms("-00", "30", "25.1")
        expected = -(0 + 30 / 60 + 25.1 / 3600)  # ≈ -0.5070°
        self.assertAlmostEqual(got, expected, delta=self.TOL)
        self.assertLess(got, 0.0)

    def test_plus_zero_is_small_positive(self):
        got = parse_declination_dms("+00", "16", "59.6")
        expected = 16 / 60 + 59.6 / 3600  # ≈ +0.2832°
        self.assertAlmostEqual(got, expected, delta=self.TOL)
        self.assertGreater(got, 0.0)

    def test_exact_zero(self):
        self.assertEqual(parse_declination_dms("+00", "00", "00.0"), 0.0)
        self.assertEqual(parse_declination_dms("00", "00", "00.0"), 0.0)
        # "-00 00 00.0" is mathematically zero. Our implementation happens
        # to return -0.0, which compares equal to 0.0 under ==; we don't
        # make any promise either way about the signed-zero bit.
        self.assertEqual(parse_declination_dms("-00", "00", "00.0"), 0.0)


# =====================================================================
# Sub-parser: parse_date_token
# =====================================================================

class TestParseDateToken(unittest.TestCase):

    def test_typical_date(self):
        self.assertEqual(parse_date_token("2035-Jan-01"), datetime(2035, 1, 1))

    def test_each_month_abbreviation(self):
        cases = [
            ("2040-Jan-15", 1), ("2040-Feb-15", 2), ("2040-Mar-15", 3),
            ("2040-Apr-15", 4), ("2040-May-15", 5), ("2040-Jun-15", 6),
            ("2040-Jul-15", 7), ("2040-Aug-15", 8), ("2040-Sep-15", 9),
            ("2040-Oct-15", 10), ("2040-Nov-15", 11), ("2040-Dec-15", 12),
        ]
        for token, expected_month in cases:
            with self.subTest(token=token):
                self.assertEqual(parse_date_token(token).month, expected_month)

    def test_leap_day(self):
        # 2036 is a leap year — a sanity check that we're not validating the
        # date ourselves; Python's datetime does it for us.
        self.assertEqual(parse_date_token("2036-Feb-29"), datetime(2036, 2, 29))

    def test_rejects_non_english_month(self):
        with self.assertRaises(ValueError):
            parse_date_token("2035-jan-01")   # lowercase — rejected
        with self.assertRaises(ValueError):
            parse_date_token("2035-JAN-01")   # uppercase — rejected
        with self.assertRaises(ValueError):
            parse_date_token("2035-Ene-01")   # Spanish — rejected

    def test_rejects_malformed(self):
        with self.assertRaises(ValueError):
            parse_date_token("2035/Jan/01")
        with self.assertRaises(ValueError):
            parse_date_token("2035-Jan")


# =====================================================================
# Integration: parse the real Horizons file
# =====================================================================

@unittest.skipUnless(
    HORIZONS_TXT.is_file(),
    f"Real Horizons file not present at {HORIZONS_TXT}",
)
class TestParseRealFile(unittest.TestCase):
    """End-to-end parse of the actual 2035-2065 ephemeris shipped in repo."""

    @classmethod
    def setUpClass(cls):
        # Silence the INFO log so the test output is clean; we verify logging
        # separately in TestParseLogging.
        logging.getLogger("coordinate_transformation.horizons_parser"
                          ).setLevel(logging.WARNING)
        cls.eph = parse(HORIZONS_TXT, log_first_entry=False)

    # ---- top-level shape ------------------------------------------------

    def test_nt_matches_day_count_2035_through_2065(self):
        # 31 years: 2035-01-01 through 2065-12-31, daily cadence.
        # Leap years in this span: 2036, 2040, 2044, 2048, 2052, 2056, 2060,
        # 2064 → 8 leap years. 31*365 + 8 = 11,323.
        self.assertEqual(self.eph.nt, 11323)

    def test_consistent_array_lengths(self):
        self.assertEqual(len(self.eph.date_strings), self.eph.nt)
        self.assertEqual(len(self.eph.datetimes), self.eph.nt)
        self.assertEqual(self.eph.subsolar_lat_deg.shape, (self.eph.nt,))

    def test_dtype_is_float64(self):
        self.assertEqual(self.eph.subsolar_lat_deg.dtype, np.float64)

    # ---- endpoints ------------------------------------------------------

    def test_first_entry(self):
        # Row 1 of the ephemeris:  2035-Jan-01 00:00 ... -23 03 57.6
        self.assertEqual(self.eph.date_strings[0], "2035-01-01")
        self.assertEqual(self.eph.datetimes[0], datetime(2035, 1, 1))
        self.assertAlmostEqual(
            self.eph.subsolar_lat_deg[0],
            -(23 + 3 / 60 + 57.6 / 3600),   # ≈ -23.0660°
            places=4,
        )

    def test_last_entry(self):
        # Last row:  2065-Dec-31 00:00 ... -23 07 48.3
        self.assertEqual(self.eph.date_strings[-1], "2065-12-31")
        self.assertEqual(self.eph.datetimes[-1], datetime(2065, 12, 31))
        self.assertAlmostEqual(
            self.eph.subsolar_lat_deg[-1],
            -(23 + 7 / 60 + 48.3 / 3600),   # ≈ -23.1301°
            places=4,
        )

    # ---- astronomically meaningful interior rows -----------------------

    def test_northern_summer_solstice_2035(self):
        # 2035-Jun-21 is the northern summer solstice — δ should be very
        # close to +ε_obliquity (≈ +23.44°). The ephemeris gives +23 25 51.7.
        idx = self.eph.date_strings.index("2035-06-21")
        self.assertAlmostEqual(
            self.eph.subsolar_lat_deg[idx],
            23 + 25 / 60 + 51.7 / 3600,     # ≈ +23.4310°
            places=4,
        )
        # Sanity: the solstice δ is within 0.02° of max obliquity.
        self.assertAlmostEqual(self.eph.subsolar_lat_deg[idx], 23.44, delta=0.02)

    def test_sign_crossing_near_march_equinox_2035(self):
        # The declination crosses zero in March. The day before crossing
        # (2035-Mar-20) is recorded as -00 30 25.1 in the file. The parser
        # must treat this as ≈ -0.5070°, not +0.5070°. This is the "-00"
        # edge case that motivates the sign convention in the module docs.
        idx = self.eph.date_strings.index("2035-03-20")
        got = self.eph.subsolar_lat_deg[idx]
        expected = -(30 / 60 + 25.1 / 3600)  # ≈ -0.5070°
        self.assertAlmostEqual(got, expected, places=4)
        self.assertLess(got, 0.0)

        # The very next day (2035-Mar-21) has δ still slightly negative:
        # the actual equinox is somewhere between Mar-21 00:00 UT and
        # Mar-22 00:00 UT in 2035.
        idx_next = self.eph.date_strings.index("2035-03-21")
        self.assertLess(self.eph.subsolar_lat_deg[idx_next], 0.0)
        # And by Mar-22 it should have flipped to positive.
        idx_after = self.eph.date_strings.index("2035-03-22")
        self.assertGreater(self.eph.subsolar_lat_deg[idx_after], 0.0)

    # ---- global physical bounds ---------------------------------------

    def test_all_values_within_obliquity(self):
        # Earth's obliquity is ~23.44°; δ can never exceed that. Allow a
        # 0.02° cushion for the decades-long nutation/drift envelope.
        lat = self.eph.subsolar_lat_deg
        self.assertTrue(np.all(np.abs(lat) <= 23.44 + 0.02),
                        msg=f"|δ| max = {np.max(np.abs(lat)):.4f}°")

    def test_mean_declination_is_small_and_slightly_positive(self):
        # Averaged over many whole years, <δ> is small but not exactly zero.
        # Earth's orbit is elliptical, with perihelion in early January
        # (northern winter) and aphelion in early July (northern summer).
        # By Kepler's 2nd law the Sun sweeps faster near perihelion, so it
        # spends fewer days below the celestial equator (~179) than above
        # (~186). That asymmetry produces a small *positive* mean δ —
        # roughly +0.3° to +0.5° depending on span.
        #
        # 31 years of daily samples should pin the mean to O(0.05°) of the
        # true long-term value; we allow a ±1° band here because the test
        # is really checking the *sign and order of magnitude* as a
        # sanity check on the parser, not a precise astronomical fit.
        mean = float(self.eph.subsolar_lat_deg.mean())
        self.assertGreater(mean, 0.0,
                           msg=f"Mean δ = {mean:+.4f}° expected positive")
        self.assertLess(abs(mean), 1.0,
                        msg=f"Mean δ = {mean:+.4f}° outside ±1° band")

    # ---- temporal continuity ------------------------------------------

    def test_datetimes_are_strictly_daily_and_sorted(self):
        dts = self.eph.datetimes
        gaps = {(dts[i + 1] - dts[i]).days for i in range(len(dts) - 1)}
        self.assertEqual(
            gaps, {1},
            msg=f"Expected exactly 1-day gaps between consecutive rows; got {gaps}",
        )

    def test_iso_date_strings_match_datetimes(self):
        for i in (0, 100, 5000, self.eph.nt - 1):
            self.assertEqual(
                self.eph.date_strings[i],
                self.eph.datetimes[i].strftime("%Y-%m-%d"),
            )


# =====================================================================
# Logging behaviour (transparency requirement)
# =====================================================================

@unittest.skipUnless(HORIZONS_TXT.is_file(), "Real Horizons file not present")
class TestParseLogging(unittest.TestCase):
    """Illeana's spec asks the parser to 'write an example of what it is
    taking'; we verify the INFO log carries the first row's date and δ."""

    def test_first_entry_logged_at_info(self):
        logger_name = "coordinate_transformation.horizons_parser"
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        handler.setLevel(logging.INFO)
        logger = logging.getLogger(logger_name)
        prev_level = logger.level
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        try:
            parse(HORIZONS_TXT, log_first_entry=True)
        finally:
            logger.removeHandler(handler)
            logger.setLevel(prev_level)

        msg = stream.getvalue()
        self.assertIn("2035-01-01", msg,
                      msg=f"First-row date missing from log: {msg!r}")
        self.assertIn("-23.066", msg,   # first 3 decimals of -23.0660
                      msg=f"First-row declination missing from log: {msg!r}")


# =====================================================================
# Error handling
# =====================================================================

class TestParseErrors(unittest.TestCase):
    """The parser should fail loudly on malformed input — silent wrong
    numbers are the worst failure mode for a scientific pipeline."""

    def _write_tmp(self, text: str) -> Path:
        # Write to a named temp file; caller is responsible for no cleanup
        # assumptions beyond the Python GC of the Path.
        tf = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, encoding="utf-8",
        )
        tf.write(text)
        tf.close()
        return Path(tf.name)

    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            parse("/no/such/horizons/file.txt")

    def test_missing_soe_marker(self):
        p = self._write_tmp(
            "header line\nmore header\n$$EOE\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse(p, log_first_entry=False)
        self.assertIn("$$SOE", str(cm.exception))

    def test_missing_eoe_marker(self):
        p = self._write_tmp(
            "header\n$$SOE\n 2035-Jan-01 00:00   18 43 02.53 -23 03 57.6\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse(p, log_first_entry=False)
        self.assertIn("$$EOE", str(cm.exception))

    def test_empty_data_section(self):
        p = self._write_tmp("header\n$$SOE\n$$EOE\n")
        with self.assertRaises(ValueError) as cm:
            parse(p, log_first_entry=False)
        self.assertIn("No ephemeris rows", str(cm.exception))

    def test_malformed_row_too_few_tokens(self):
        p = self._write_tmp(
            "header\n$$SOE\n 2035-Jan-01 00:00   18 43 02.53\n$$EOE\n"
        )
        with self.assertRaises(ValueError) as cm:
            parse(p, log_first_entry=False)
        # Error should point at the offending line and say what was expected.
        self.assertIn("8", str(cm.exception))  # "at least 8 tokens" message

    def test_minimal_synthetic_roundtrip(self):
        # Construct a tiny three-day fixture and verify we parse it back
        # exactly — a confidence check on the parser independent of the
        # real Horizons file, useful if that file ever moves.
        p = self._write_tmp(
            "irrelevant header\n"
            "more header lines\n"
            "$$SOE\n"
            " 2040-Jan-01 00:00     18 43 02.53 -23 03 57.6\n"
            " 2040-Jan-02 00:00     18 47 27.49 +00 00 00.0\n"
            " 2040-Jan-03 00:00     18 51 52.15 +12 34 56.7\n"
            "$$EOE\n"
            "trailing footer\n"
        )
        eph = parse(p, log_first_entry=False)
        self.assertEqual(eph.nt, 3)
        self.assertEqual(
            eph.date_strings, ["2040-01-01", "2040-01-02", "2040-01-03"],
        )
        self.assertAlmostEqual(
            float(eph.subsolar_lat_deg[0]),
            -(23 + 3 / 60 + 57.6 / 3600), places=4,
        )
        self.assertEqual(float(eph.subsolar_lat_deg[1]), 0.0)
        self.assertAlmostEqual(
            float(eph.subsolar_lat_deg[2]),
            12 + 34 / 60 + 56.7 / 3600, places=4,
        )
        self.assertIsInstance(eph, HorizonsEphemeris)


if __name__ == "__main__":
    unittest.main(verbosity=2)
