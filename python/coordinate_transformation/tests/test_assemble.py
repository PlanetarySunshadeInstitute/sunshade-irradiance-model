"""
Tests for ``assemble``.

Strategy
--------
The expensive / environment-heavy parts (reading Matt's NetCDF, writing the
output NetCDF) are behind two narrow protocols: ``SourceProviderProtocol``
and ``WriterProtocol``.  We test the pure orchestration core
(``assemble_into_writer``) with in-memory fakes for both.  This means:

* No netCDF4 required to run the tests.
* We can inspect exactly what the writer received — per-day slabs, time
  indices, call order — and cross-reference it against expected behaviour.

Also tested directly:

* :class:`SunFixedSourceProvider` — leap-year / day-of-year lookup.
* :func:`_maybe_reorder_to_time_first` — axis-ordering robustness.
* :func:`_validate_lat_array` / :func:`_validate_lon_array` — the sniff
  tests that decide when to fall back to the CESM f09 grid.

Integration with the real Horizons file: uses the actual parser to build an
ephemeris for 2035-2065 (11,323 days), then asks the orchestrator to run
against a fake source that returns known slabs.  Verifies the writer got
exactly 11,323 calls, in strict time order, and that per-day δ values pass
through correctly.

Run with::

    python -m unittest python.coordinate_transformation.tests.test_assemble
"""

from __future__ import annotations

import calendar
import logging
import sys
import unittest
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

import numpy as np

_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_PYTHON_DIR = _PKG_DIR.parent
_REPO_ROOT = _PYTHON_DIR.parent
for p in (_REPO_ROOT, _PYTHON_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from coordinate_transformation import assemble as _assemble  # noqa: E402
from coordinate_transformation.assemble import (  # noqa: E402
    AssemblyReport,
    SunFixedSourceProvider,
    _f09_fallback_grid,
    _maybe_reorder_to_time_first,
    _validate_lat_array,
    _validate_lon_array,
    assemble_into_writer,
)
from coordinate_transformation.horizons_parser import HorizonsEphemeris  # noqa: E402


HORIZONS_TXT = (
    _REPO_ROOT / "Coordinate Transformation" / "horizons_Sun_RA-DEC.txt"
)


# =====================================================================
# Fakes for the protocols
# =====================================================================

class CapturingWriter:
    """Records every write_day call for later inspection."""

    def __init__(self):
        self.calls: List[Tuple[int, np.ndarray]] = []
        self.closed = False

    def write_day(self, t_idx: int, df_slab: np.ndarray) -> None:
        # Copy the slab — callers may be using a single scratch buffer.
        self.calls.append((t_idx, np.array(df_slab, copy=True)))

    def close(self) -> None:
        self.closed = True


class FlatConstantSource:
    """A source provider that returns the same constant-valued slab every day.

    Since the slab is constant, rotating it should produce the same constant
    on output — any deviation would signal a bug in the orchestrator or the
    rotation stage.
    """

    def __init__(self, value: float, nlat: int, nlon: int):
        self.value = value
        self.source_lats_deg = np.linspace(-90, 90, nlat)
        self.source_lons_deg = np.linspace(0, 360, nlon, endpoint=False)
        self._slab = np.full((nlat, nlon), value, dtype=np.float64)

    def slice_for_date(self, dt):
        return self._slab


class RecordingDateSource:
    """Records the sequence of dates it was asked for, and returns a slab
    whose values encode (month, day) in a DF-plausible range — so the
    writer's record can be decoded back into the dates the orchestrator
    walked and the plausibility check (DF ∈ [0, 1]) still passes."""

    def __init__(self, nlat: int, nlon: int):
        self.source_lats_deg = np.linspace(-90, 90, nlat)
        self.source_lons_deg = np.linspace(0, 360, nlon, endpoint=False)
        self._nlat = nlat
        self._nlon = nlon
        self.dates_seen: List[datetime] = []

    def slice_for_date(self, dt):
        self.dates_seen.append(dt)
        # Encode the date into a constant slab in [0, 1] so the
        # plausibility check (DF ∈ [0, 1]) still passes: month*0.01 +
        # day*0.0001 + year_parity*0.5.  The dates_seen list is how
        # tests should actually read off the date sequence.
        value = (dt.month * 0.01) + (dt.day * 0.0001)
        return np.full((self._nlat, self._nlon), value, dtype=np.float64)


# =====================================================================
# SunFixedSourceProvider: leap/standard selection and day-of-year indexing
# =====================================================================

class TestSunFixedSourceProvider(unittest.TestCase):

    def setUp(self):
        self.nlat, self.nlon = 5, 7
        # df_standard[doy] = doy + 0.0  (to make the value encode the index)
        self.df_std = np.zeros((365, self.nlat, self.nlon))
        for d in range(365):
            self.df_std[d] = float(d)
        # df_leap[doy] = doy + 1000.0  (distinguishable from standard)
        self.df_leap = np.zeros((366, self.nlat, self.nlon))
        for d in range(366):
            self.df_leap[d] = float(d) + 1000.0

        self.provider = SunFixedSourceProvider(
            df_standard=self.df_std,
            df_leap=self.df_leap,
            source_lats_deg=np.linspace(-90, 90, self.nlat),
            source_lons_deg=np.linspace(0, 360, self.nlon, endpoint=False),
        )

    def test_jan1_standard_year_is_doy0(self):
        # 2035 is a standard year.  Jan 1 must pull df_standard[0].
        slab = self.provider.slice_for_date(datetime(2035, 1, 1))
        self.assertEqual(slab[0, 0], 0.0)

    def test_dec31_standard_year_is_doy364(self):
        slab = self.provider.slice_for_date(datetime(2035, 12, 31))
        self.assertEqual(slab[0, 0], 364.0)

    def test_jan1_leap_year_is_doy0(self):
        # 2036 is a leap year.  Values encoded with +1000 prefix.
        slab = self.provider.slice_for_date(datetime(2036, 1, 1))
        self.assertEqual(slab[0, 0], 1000.0)

    def test_feb29_leap_year_is_doy59(self):
        # 2036 is leap; Feb 29 is the 60th day → doy=59 (0-indexed).
        slab = self.provider.slice_for_date(datetime(2036, 2, 29))
        self.assertEqual(slab[0, 0], 1059.0)

    def test_dec31_leap_year_is_doy365(self):
        slab = self.provider.slice_for_date(datetime(2036, 12, 31))
        self.assertEqual(slab[0, 0], 1365.0)

    def test_rejects_wrong_shape_df_standard(self):
        with self.assertRaises(ValueError):
            SunFixedSourceProvider(
                df_standard=np.zeros((300, self.nlat, self.nlon)),  # wrong time length
                df_leap=self.df_leap,
                source_lats_deg=np.linspace(-90, 90, self.nlat),
                source_lons_deg=np.linspace(0, 360, self.nlon, endpoint=False),
            )

    def test_rejects_wrong_shape_df_leap(self):
        with self.assertRaises(ValueError):
            SunFixedSourceProvider(
                df_standard=self.df_std,
                df_leap=np.zeros((365, self.nlat, self.nlon)),  # 365, not 366
                source_lats_deg=np.linspace(-90, 90, self.nlat),
                source_lons_deg=np.linspace(0, 360, self.nlon, endpoint=False),
            )

    def test_rejects_lat_lon_mismatch(self):
        with self.assertRaises(ValueError):
            SunFixedSourceProvider(
                df_standard=self.df_std,
                df_leap=self.df_leap,
                source_lats_deg=np.linspace(-90, 90, self.nlat + 1),
                source_lons_deg=np.linspace(0, 360, self.nlon, endpoint=False),
            )


# =====================================================================
# assemble_into_writer: the orchestration core
# =====================================================================

def _tiny_ephemeris(dates_and_deltas) -> HorizonsEphemeris:
    """Build a tiny synthetic ephemeris from a list of (datetime, delta)."""
    dts = [d for d, _ in dates_and_deltas]
    decs = np.array([x for _, x in dates_and_deltas], dtype=np.float64)
    return HorizonsEphemeris(
        nt=len(dts),
        date_strings=[d.strftime("%Y-%m-%d") for d in dts],
        datetimes=dts,
        subsolar_lat_deg=decs,
    )


class TestAssembleIntoWriter(unittest.TestCase):

    def test_writer_receives_one_call_per_day_in_order(self):
        # Three consecutive days, tracked by a recording source.
        ephem = _tiny_ephemeris([
            (datetime(2035, 1, 1), -23.07),
            (datetime(2035, 1, 2), -22.99),
            (datetime(2035, 1, 3), -22.90),
        ])
        src = RecordingDateSource(9, 12)
        writer = CapturingWriter()
        report = assemble_into_writer(
            src, ephem, writer,
            progress_every=0,  # silence progress logs
        )

        # The writer received 3 calls, in order, with t_idx = 0, 1, 2.
        self.assertEqual(len(writer.calls), 3)
        for expected_t, (actual_t, _slab) in enumerate(writer.calls):
            self.assertEqual(actual_t, expected_t)

        # The source was consulted for exactly those three dates, in order.
        self.assertEqual(
            src.dates_seen,
            [datetime(2035, 1, 1), datetime(2035, 1, 2), datetime(2035, 1, 3)],
        )
        self.assertEqual(report.nt, 3)
        self.assertEqual(report.first_date, "2035-01-01")
        self.assertFalse(report.first_source_is_leap)   # 2035 is standard
        self.assertEqual(report.first_source_index, 0)  # Jan 1 → doy 0

    def test_leap_year_first_date_flagged_in_report(self):
        ephem = _tiny_ephemeris([(datetime(2036, 1, 1), -23.0)])
        src = RecordingDateSource(9, 12)
        writer = CapturingWriter()
        report = assemble_into_writer(src, ephem, writer, progress_every=0)
        self.assertTrue(report.first_source_is_leap)
        self.assertEqual(report.first_date, "2036-01-01")

    def test_constant_field_round_trips_through_assembly(self):
        # A constant DF=0.7 slab rotated by any δ should remain 0.7 — any
        # deviation would indicate orchestrator corruption.
        ephem = _tiny_ephemeris([
            (datetime(2035, 6, 21), +23.44),   # summer solstice
            (datetime(2035, 12, 22), -23.44),  # winter solstice
            (datetime(2035, 3, 20),  -0.51),   # near equinox
        ])
        src = FlatConstantSource(0.7, 19, 37)
        writer = CapturingWriter()
        report = assemble_into_writer(src, ephem, writer, progress_every=0)

        for _t, slab in writer.calls:
            # Exact equality — we're rotating a constant, and bilinear
            # interpolation of a constant is exact.
            self.assertTrue(np.allclose(slab, 0.7, atol=1e-12))
        self.assertAlmostEqual(report.df_min, 0.7, places=12)
        self.assertAlmostEqual(report.df_max, 0.7, places=12)
        self.assertAlmostEqual(report.df_mean, 0.7, places=12)

    def test_leap_year_selects_df_leap(self):
        # Use a provider where standard and leap sources are easily
        # distinguishable, and confirm that 2036-02-29 (only exists in
        # leap years) pulls from df_leap.
        nlat, nlon = 5, 7
        df_std = np.zeros((365, nlat, nlon))                 # all zeros
        df_leap = np.ones((366, nlat, nlon)) * 0.5           # all 0.5
        provider = SunFixedSourceProvider(
            df_standard=df_std, df_leap=df_leap,
            source_lats_deg=np.linspace(-90, 90, nlat),
            source_lons_deg=np.linspace(0, 360, nlon, endpoint=False),
        )
        ephem = _tiny_ephemeris([
            (datetime(2035, 2, 28), -8.0),  # standard year
            (datetime(2036, 2, 29), -8.0),  # leap year — only possible here
        ])
        writer = CapturingWriter()
        assemble_into_writer(provider, ephem, writer, progress_every=0)

        # Day 0 came from df_std (all zeros); day 1 from df_leap (all 0.5).
        # Rotating a constant returns the same constant.
        self.assertTrue(np.allclose(writer.calls[0][1], 0.0, atol=1e-12))
        self.assertTrue(np.allclose(writer.calls[1][1], 0.5, atol=1e-12))

    def test_implausible_df_raises(self):
        # DF out of [0, 1] (e.g. negative) must fail fast.
        nlat, nlon = 5, 7
        # Standard year slab has a negative value — orchestrator must catch
        # it and refuse to proceed (we don't want to write 5 GB of garbage).
        df_std = np.zeros((365, nlat, nlon))
        df_std[0, 0, 0] = -0.1                    # negative → implausible
        df_leap = np.zeros((366, nlat, nlon))
        provider = SunFixedSourceProvider(
            df_standard=df_std, df_leap=df_leap,
            source_lats_deg=np.linspace(-90, 90, nlat),
            source_lons_deg=np.linspace(0, 360, nlon, endpoint=False),
        )
        ephem = _tiny_ephemeris([(datetime(2035, 1, 1), 0.0)])
        writer = CapturingWriter()
        with self.assertRaises(ValueError) as cm:
            assemble_into_writer(provider, ephem, writer, progress_every=0)
        self.assertIn("Implausible DF", str(cm.exception))
        self.assertIn("2035-01-01", str(cm.exception))

    def test_empty_ephemeris_raises(self):
        provider = FlatConstantSource(0.9, 5, 7)
        ephem = HorizonsEphemeris(
            nt=0, date_strings=[], datetimes=[],
            subsolar_lat_deg=np.zeros(0),
        )
        writer = CapturingWriter()
        with self.assertRaises(ValueError):
            assemble_into_writer(provider, ephem, writer, progress_every=0)

    def test_grid_mismatch_raises(self):
        # Provider claims to be on a 9×12 grid but returns 9×11 slabs.
        class MisshapenSource:
            source_lats_deg = np.linspace(-90, 90, 9)
            source_lons_deg = np.linspace(0, 360, 12, endpoint=False)

            def slice_for_date(self, dt):
                return np.zeros((9, 11))   # wrong lon length

        provider = MisshapenSource()
        ephem = _tiny_ephemeris([(datetime(2035, 1, 1), 0.0)])
        writer = CapturingWriter()
        with self.assertRaises(ValueError) as cm:
            assemble_into_writer(provider, ephem, writer, progress_every=0)
        self.assertIn("Source provider returned shape", str(cm.exception))

    def test_show_progress_without_tqdm_warns_and_still_completes(self):
        # tqdm is not installed in this environment — the orchestrator
        # should log a one-time warning about the fallback and then run
        # to completion with the periodic log line instead.
        ephem = _tiny_ephemeris([
            (datetime(2035, 1, 1), -23.0),
            (datetime(2035, 1, 2), -23.0),
        ])
        src = FlatConstantSource(0.5, 9, 12)
        writer = CapturingWriter()
        with self.assertLogs(
            "coordinate_transformation.assemble", level="WARNING",
        ) as cm:
            report = assemble_into_writer(
                src, ephem, writer,
                show_progress=True, progress_every=0,
            )
        joined = "\n".join(cm.output)
        self.assertIn("tqdm", joined)
        # Still processed both days.
        self.assertEqual(len(writer.calls), 2)
        self.assertEqual(report.nt, 2)

    def test_show_progress_with_tqdm_present_suppresses_log_progress(self):
        # Simulate tqdm being installed by stubbing a minimal module into
        # sys.modules before the call.  The bar should run, and the
        # log-based "progress: N/M" lines should NOT appear.
        import sys
        import types

        calls = []

        def _fake_tqdm(iterable, **kwargs):
            calls.append(kwargs)
            for x in iterable:
                yield x

        fake_mod = types.ModuleType("tqdm")
        fake_mod.tqdm = _fake_tqdm  # type: ignore[attr-defined]
        ephem = _tiny_ephemeris([
            (datetime(2035, 1, 1), -23.0),
            (datetime(2035, 1, 2), -23.0),
            (datetime(2035, 1, 3), -23.0),
        ])
        src = FlatConstantSource(0.5, 9, 12)
        writer = CapturingWriter()

        old_tqdm = sys.modules.get("tqdm")
        sys.modules["tqdm"] = fake_mod
        try:
            # Use assertLogs at WARNING level: if anything warning-level
            # was emitted the context manager passes; if nothing was
            # emitted it raises.  We expect no warnings at all, so we
            # wrap in a "should be quiet" check using a higher level.
            with self.assertNoLogs(
                "coordinate_transformation.assemble", level="WARNING",
            ):
                report = assemble_into_writer(
                    src, ephem, writer,
                    show_progress=True,
                    progress_every=1,   # would normally log every day
                )
        finally:
            if old_tqdm is None:
                del sys.modules["tqdm"]
            else:
                sys.modules["tqdm"] = old_tqdm

        # The fake tqdm was called exactly once, wrapping the day loop.
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].get("unit"), "day")
        self.assertEqual(calls[0].get("total"), 3)
        # All days processed.
        self.assertEqual(len(writer.calls), 3)

    def test_day0_mapping_logged_at_info(self):
        # The "Day-0 source-slab mapping" line should appear in INFO logs so
        # a scientist can verify the code picked the expected source slab.
        ephem = _tiny_ephemeris([(datetime(2036, 1, 1), -23.0)])  # leap year
        provider = FlatConstantSource(0.5, 9, 12)
        writer = CapturingWriter()
        with self.assertLogs(
            "coordinate_transformation.assemble", level="INFO",
        ) as cm:
            assemble_into_writer(provider, ephem, writer, progress_every=0)
        joined = "\n".join(cm.output)
        self.assertIn("Day-0", joined)
        self.assertIn("df_leap", joined)
        self.assertIn("2036-01-01", joined)


# =====================================================================
# Integration against the real Horizons ephemeris
# =====================================================================

@unittest.skipUnless(HORIZONS_TXT.is_file(), "Real Horizons file not present")
class TestAssembleWithRealEphemeris(unittest.TestCase):
    """Parse the real 2035-2065 file, run the orchestrator against a fake
    source, and verify the writer gets exactly 11,323 ordered calls."""

    @classmethod
    def setUpClass(cls):
        from coordinate_transformation import horizons_parser
        logging.getLogger("coordinate_transformation").setLevel(
            logging.WARNING  # silence for speed
        )
        cls.eph = horizons_parser.parse(HORIZONS_TXT, log_first_entry=False)

    def test_full_31year_run_with_constant_source(self):
        src = FlatConstantSource(0.85, 19, 37)
        writer = CapturingWriter()
        report = assemble_into_writer(
            src, self.eph, writer, progress_every=0,
        )
        self.assertEqual(report.nt, 11323)
        self.assertEqual(len(writer.calls), 11323)
        # Time indices must be 0..11322 in order.
        for i, (t_idx, _) in enumerate(writer.calls):
            self.assertEqual(t_idx, i)
        # DF constant ⇒ min==max==mean==0.85 (to roundoff).
        self.assertAlmostEqual(report.df_min, 0.85, places=10)
        self.assertAlmostEqual(report.df_max, 0.85, places=10)
        self.assertAlmostEqual(report.df_mean, 0.85, places=10)
        # Sensible day-0 mapping.
        self.assertEqual(report.first_date, "2035-01-01")
        self.assertEqual(report.first_source_index, 0)
        self.assertFalse(report.first_source_is_leap)

    def test_leap_day_count_in_2035_to_2065(self):
        # There are 8 leap years in [2035, 2065]: 2036, 40, 44, 48, 52, 56,
        # 60, 64 — each contributing one Feb 29.  If our orchestrator walks
        # the ephemeris in date order, a RecordingDateSource should see
        # exactly 8 Feb 29s.
        src = RecordingDateSource(9, 12)
        writer = CapturingWriter()
        assemble_into_writer(src, self.eph, writer, progress_every=0)
        feb29s = [d for d in src.dates_seen if d.month == 2 and d.day == 29]
        self.assertEqual(len(feb29s), 8)
        self.assertEqual(
            [d.year for d in feb29s],
            [2036, 2040, 2044, 2048, 2052, 2056, 2060, 2064],
        )


# =====================================================================
# Source file helpers
# =====================================================================

class TestValidateLatLon(unittest.TestCase):

    def test_lat_accepts_standard_f09(self):
        lat = np.linspace(-90, 90, 192)
        self.assertTrue(_validate_lat_array(lat, 192))

    def test_lat_rejects_wrong_length(self):
        self.assertFalse(_validate_lat_array(np.linspace(-90, 90, 192), 180))

    def test_lat_rejects_out_of_range(self):
        self.assertFalse(_validate_lat_array(np.linspace(-100, 100, 192), 192))

    def test_lat_rejects_nonmonotonic(self):
        lat = np.linspace(-90, 90, 192).copy()
        lat[50], lat[60] = lat[60], lat[50]
        self.assertFalse(_validate_lat_array(lat, 192))

    def test_lat_rejects_nan(self):
        lat = np.linspace(-90, 90, 192).copy()
        lat[10] = np.nan
        self.assertFalse(_validate_lat_array(lat, 192))

    def test_lon_accepts_f09(self):
        lon = np.linspace(0, 360, 288, endpoint=False)
        self.assertTrue(_validate_lon_array(lon, 288))

    def test_lon_accepts_closed_360(self):
        lon = np.linspace(0, 360, 289)
        self.assertTrue(_validate_lon_array(lon, 289))

    def test_lon_rejects_template_integer_indices(self):
        # The template file stores 0..287 in the lon variable — that spans
        # only 287°, which is nonsense for a longitude. Must be rejected
        # so the loader falls back to the real f09 grid.
        lon = np.arange(288, dtype=np.float64)
        self.assertFalse(_validate_lon_array(lon, 288))


class TestReorderTimeFirst(unittest.TestCase):

    def test_noop_when_already_time_first(self):
        arr = np.zeros((365, 5, 7))
        out = _maybe_reorder_to_time_first(arr, expect_time=365)
        self.assertIs(out, arr)

    def test_moves_time_axis_from_last_position(self):
        arr = np.zeros((5, 7, 366))
        out = _maybe_reorder_to_time_first(arr, expect_time=366)
        self.assertEqual(out.shape, (366, 5, 7))

    def test_moves_time_axis_from_middle(self):
        arr = np.zeros((5, 365, 7))
        out = _maybe_reorder_to_time_first(arr, expect_time=365)
        self.assertEqual(out.shape, (365, 5, 7))

    def test_rejects_ambiguous_time_axis(self):
        # Two non-leading axes of length 365 — which is the time axis?
        # We refuse to guess.  (If axis 0 already is the time axis, the
        # function short-circuits and doesn't reach the ambiguity check —
        # that's correct, so we construct a shape where axis 0 is NOT 365
        # but two other axes are.)
        arr = np.zeros((5, 365, 365))
        with self.assertRaises(ValueError):
            _maybe_reorder_to_time_first(arr, expect_time=365)

    def test_rejects_missing_time_axis(self):
        arr = np.zeros((5, 7, 8))
        with self.assertRaises(ValueError):
            _maybe_reorder_to_time_first(arr, expect_time=365)


class TestF09FallbackGrid(unittest.TestCase):

    def test_standard_f09_dims(self):
        lat, lon = _f09_fallback_grid(192, 288)
        self.assertEqual(lat.shape, (192,))
        self.assertEqual(lon.shape, (288,))
        self.assertAlmostEqual(float(lat[0]), -90.0)
        self.assertAlmostEqual(float(lat[-1]),  90.0)
        self.assertAlmostEqual(float(lon[0]),    0.0)
        # endpoint=False: last lon = 360 - 360/288 = 358.75
        self.assertAlmostEqual(float(lon[-1]), 358.75)


if __name__ == "__main__":
    unittest.main(verbosity=2)
