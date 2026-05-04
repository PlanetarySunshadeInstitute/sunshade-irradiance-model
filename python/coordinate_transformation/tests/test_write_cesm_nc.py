"""
Tests for ``write_cesm_nc``.

Two layers:

1. **Pure-helper tests** — exercise ``encode_cf_time``,
   ``build_global_attrs``, ``resolve_atomic_paths``, and the
   ``WriterOptions`` dataclass.  These run in any environment without
   netCDF4 installed.

2. **End-to-end integration tests** — actually open a tiny output file
   with the writer, write a few days of synthetic DF, close, then
   re-open and verify dimensions, variables, attributes, dtypes, and
   data values.  These ``unittest.skipUnless`` netCDF4 is importable.

Run with::

    python -m unittest python.coordinate_transformation.tests.test_write_cesm_nc
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_PYTHON_DIR = _PKG_DIR.parent
_REPO_ROOT = _PYTHON_DIR.parent
for p in (_REPO_ROOT, _PYTHON_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from coordinate_transformation.horizons_parser import HorizonsEphemeris  # noqa: E402
from coordinate_transformation.write_cesm_nc import (  # noqa: E402
    CF_CONVENTIONS,
    DF_FILL_VALUE,
    DIM_LAT,
    DIM_LON,
    DIM_TIME,
    STRING_LEN,
    VAR_DATE,
    VAR_DF,
    VAR_LAT,
    VAR_LON,
    VAR_TIME,
    WriterOptions,
    build_global_attrs,
    encode_cf_time,
    resolve_atomic_paths,
)


def _tiny_ephemeris(n: int = 5, start: datetime = datetime(2035, 1, 1)) -> HorizonsEphemeris:
    """Build a synthetic ephemeris of ``n`` consecutive daily entries."""
    dts = [start + timedelta(days=i) for i in range(n)]
    return HorizonsEphemeris(
        nt=n,
        date_strings=[d.strftime("%Y-%m-%d") for d in dts],
        datetimes=dts,
        # Realistic-ish daily declinations near January.
        subsolar_lat_deg=np.linspace(-23.0, -22.5, n).astype(np.float64),
    )


# =====================================================================
# encode_cf_time
# =====================================================================

class TestEncodeCfTime(unittest.TestCase):

    def test_default_reference_is_first_day_at_midnight(self):
        dts = [datetime(2035, 1, 1), datetime(2035, 1, 2), datetime(2035, 1, 3)]
        values, units = encode_cf_time(dts)
        np.testing.assert_array_equal(values, [0.0, 1.0, 2.0])
        self.assertEqual(units, "days since 2035-01-01 00:00:00")

    def test_explicit_reference(self):
        dts = [datetime(2035, 6, 21), datetime(2035, 6, 22)]
        ref = datetime(2035, 1, 1)
        values, units = encode_cf_time(dts, reference=ref)
        # Jan 1 → Jun 21 in 2035 (standard year): 31+28+31+30+31+20 = 171 days
        self.assertAlmostEqual(float(values[0]), 171.0)
        self.assertAlmostEqual(float(values[1]), 172.0)
        self.assertEqual(units, "days since 2035-01-01 00:00:00")

    def test_handles_leap_year_correctly(self):
        # 2036 is a leap year.  Jan 1 2036 → Jan 1 2037 should be 366 days.
        ref = datetime(2036, 1, 1)
        values, _ = encode_cf_time([datetime(2037, 1, 1)], reference=ref)
        self.assertAlmostEqual(float(values[0]), 366.0)

    def test_subdaily_precision_preserved(self):
        ref = datetime(2035, 1, 1)
        dts = [datetime(2035, 1, 1, 12, 0, 0)]      # noon, half a day
        values, _ = encode_cf_time(dts, reference=ref)
        self.assertAlmostEqual(float(values[0]), 0.5, places=10)

    def test_full_31_year_span_endpoints(self):
        # Sanity check on the actual production case.
        dts = [datetime(2035, 1, 1), datetime(2065, 12, 31)]
        values, units = encode_cf_time(dts)
        # 31 years includes leap days 2036, 40, 44, 48, 52, 56, 60, 64 → 8.
        # Standard span: 2065-12-31 minus 2035-01-01 = 31 years - 1 day, but
        # exact day count is what matters; let's just check it's positive
        # and matches Python's own arithmetic.
        expected_days = (datetime(2065, 12, 31) - datetime(2035, 1, 1)).days
        self.assertEqual(int(values[1]), expected_days)
        self.assertEqual(units, "days since 2035-01-01 00:00:00")

    def test_rejects_empty_input(self):
        with self.assertRaises(ValueError):
            encode_cf_time([])


# =====================================================================
# build_global_attrs
# =====================================================================

class TestBuildGlobalAttrs(unittest.TestCase):

    def test_required_keys_present(self):
        eph = _tiny_ephemeris()
        attrs = build_global_attrs(eph)
        for key in (
            "Conventions", "title", "institution",
            "source", "references", "history", "comment",
        ):
            self.assertIn(key, attrs)

    def test_conventions_is_cf18(self):
        attrs = build_global_attrs(_tiny_ephemeris())
        self.assertEqual(attrs["Conventions"], CF_CONVENTIONS)
        self.assertTrue(attrs["Conventions"].startswith("CF-"))

    def test_history_includes_period_and_paths(self):
        eph = _tiny_ephemeris()
        attrs = build_global_attrs(
            eph,
            source_nc_path="/some/path/source.nc",
            horizons_path="/some/path/horizons.txt",
        )
        history = attrs["history"]
        self.assertIn("source NetCDF", history)
        self.assertIn("/some/path/source.nc", history)
        self.assertIn("ephemeris", history)
        self.assertIn("/some/path/horizons.txt", history)
        self.assertIn(eph.date_strings[0], history)
        self.assertIn(eph.date_strings[-1], history)

    def test_history_works_without_optional_paths(self):
        # Optional paths omitted: history still records the period.
        attrs = build_global_attrs(_tiny_ephemeris())
        self.assertIn("daily snapshots", attrs["history"])

    def test_extra_history_appended(self):
        attrs = build_global_attrs(
            _tiny_ephemeris(), extra_history="custom postscript here",
        )
        self.assertIn("custom postscript here", attrs["history"])


# =====================================================================
# resolve_atomic_paths
# =====================================================================

class TestResolveAtomicPaths(unittest.TestCase):

    def test_partial_suffix_appended_to_filename(self):
        final, partial = resolve_atomic_paths("/tmp/output.nc")
        self.assertEqual(str(final), "/tmp/output.nc")
        self.assertEqual(str(partial), "/tmp/output.nc.partial")

    def test_works_with_path_object(self):
        final, partial = resolve_atomic_paths(Path("/tmp/output.nc"))
        self.assertIsInstance(final, Path)
        self.assertIsInstance(partial, Path)
        self.assertEqual(partial.name, "output.nc.partial")

    def test_partial_is_in_same_directory(self):
        final, partial = resolve_atomic_paths("/some/dir/out.nc")
        self.assertEqual(final.parent, partial.parent)


# =====================================================================
# WriterOptions
# =====================================================================

class TestWriterOptions(unittest.TestCase):

    def test_defaults(self):
        opt = WriterOptions()
        self.assertTrue(opt.compress)
        self.assertEqual(opt.compress_level, 4)
        self.assertEqual(opt.nc_format, "NETCDF4")
        self.assertTrue(opt.remove_partial_on_error)


# =====================================================================
# End-to-end integration — requires netCDF4
# =====================================================================

try:
    import netCDF4 as _nc4   # noqa: F401
    HAS_NETCDF4 = True
except ImportError:
    HAS_NETCDF4 = False


@unittest.skipUnless(HAS_NETCDF4, "netCDF4 not installed in this environment")
class TestCesmNcWriterIntegration(unittest.TestCase):
    """Open the writer, write a few days, close, re-open, and verify
    everything that ``ncdump -h`` would show plus the actual data."""

    def setUp(self):
        from coordinate_transformation.write_cesm_nc import CesmNcWriter
        self.WriterCls = CesmNcWriter

        # Tiny grid, 5 days
        self.nlat, self.nlon, self.nt = 9, 12, 5
        self.lats = np.linspace(-90, 90, self.nlat)
        self.lons = np.linspace(0, 360, self.nlon, endpoint=False)
        self.eph = _tiny_ephemeris(self.nt)

        self.tmpdir = tempfile.mkdtemp(prefix="cesm_nc_test_")
        self.out_path = Path(self.tmpdir) / "out.nc"

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _write_simple_file(self):
        """Write a small file with df_slab[t, i, j] = t * 0.1 + 0.01 — a
        per-day-distinguishable but plausible (DF ∈ [0,1]) value."""
        with self.WriterCls(
            path=self.out_path,
            source_lats_deg=self.lats,
            source_lons_deg=self.lons,
            ephemeris=self.eph,
        ) as w:
            for t in range(self.nt):
                slab = np.full((self.nlat, self.nlon), t * 0.1 + 0.01)
                w.write_day(t, slab)
        return self.out_path

    # ---- Filesystem behaviour ------------------------------------------

    def test_partial_renamed_to_final_on_clean_close(self):
        path = self._write_simple_file()
        self.assertTrue(path.is_file())
        self.assertFalse(path.with_name(path.name + ".partial").exists())

    def test_partial_removed_on_exception_in_context(self):
        # Force an exception inside the writer context — partial file
        # must be cleaned up.
        with self.assertRaises(RuntimeError):
            with self.WriterCls(
                path=self.out_path,
                source_lats_deg=self.lats,
                source_lons_deg=self.lons,
                ephemeris=self.eph,
            ) as w:
                w.write_day(0, np.full((self.nlat, self.nlon), 0.5))
                raise RuntimeError("boom")
        self.assertFalse(self.out_path.exists())
        self.assertFalse(self.out_path.with_name(
            self.out_path.name + ".partial").exists())

    # ---- Schema ---------------------------------------------------------

    def test_dimensions_written(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            self.assertEqual(len(ds.dimensions[DIM_TIME]), self.nt)
            self.assertEqual(len(ds.dimensions[DIM_LAT]), self.nlat)
            self.assertEqual(len(ds.dimensions[DIM_LON]), self.nlon)

    def test_coordinate_variables_written(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            np.testing.assert_array_equal(ds[VAR_LAT][:], self.lats)
            np.testing.assert_array_equal(ds[VAR_LON][:], self.lons)
            # time: 0..nt-1 (days since first day)
            np.testing.assert_array_equal(
                ds[VAR_TIME][:], np.arange(self.nt, dtype=np.float64),
            )

    def test_cf_time_attributes(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            t = ds[VAR_TIME]
            # Option A: DF slabs are labelled at noon UT; reference epoch
            # is noon of the first day so values stay integer-valued days.
            self.assertEqual(t.units, "days since 2035-01-01 12:00:00")
            self.assertEqual(t.calendar, "gregorian")
            self.assertEqual(t.standard_name, "time")
            self.assertEqual(t.axis, "T")
            self.assertIn("12:00 UT", t.long_name)

    def test_cf_lat_lon_attributes(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            self.assertEqual(ds[VAR_LAT].units, "degrees_north")
            self.assertEqual(ds[VAR_LAT].standard_name, "latitude")
            self.assertEqual(ds[VAR_LAT].axis, "Y")
            self.assertEqual(ds[VAR_LON].units, "degrees_east")
            self.assertEqual(ds[VAR_LON].standard_name, "longitude")
            self.assertEqual(ds[VAR_LON].axis, "X")

    def test_date_strings_round_trip(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            chars = ds[VAR_DATE][:]   # (nt, STRING_LEN) char array
            decoded = _nc4.chartostring(chars)
            np.testing.assert_array_equal(
                np.asarray(decoded, dtype="U10"),
                np.asarray(self.eph.date_strings, dtype="U10"),
            )

    def test_df_values_round_trip(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            df = ds[VAR_DF][:]
            self.assertEqual(df.shape, (self.nt, self.nlat, self.nlon))
            for t in range(self.nt):
                expected = t * 0.1 + 0.01
                np.testing.assert_allclose(df[t], expected, atol=1e-12)

    def test_df_attributes(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            df = ds[VAR_DF]
            self.assertEqual(df.long_name, "dimming factor")
            self.assertEqual(df.units, "1")
            np.testing.assert_array_equal(df.valid_range, [0.0, 1.0])
            self.assertAlmostEqual(float(df._FillValue), DF_FILL_VALUE)

    def test_global_attributes(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            self.assertEqual(ds.Conventions, CF_CONVENTIONS)
            self.assertIn("Planetary Sunshade", ds.title)
            self.assertEqual(ds.institution, "Planetary Sunshade Institute")
            self.assertIn("created by", ds.history)
            self.assertIn(self.eph.date_strings[0], ds.history)
            self.assertIn(self.eph.date_strings[-1], ds.history)

    def test_df_chunking(self):
        path = self._write_simple_file()
        with _nc4.Dataset(str(path), "r") as ds:
            df = ds[VAR_DF]
            # Per-day chunking: first chunk dim is 1, others are full axis.
            chunking = df.chunking()
            self.assertEqual(chunking, [1, self.nlat, self.nlon])

    # ---- Validation ----------------------------------------------------

    def test_rejects_writes_outside_t_range(self):
        with self.WriterCls(
            path=self.out_path,
            source_lats_deg=self.lats,
            source_lons_deg=self.lons,
            ephemeris=self.eph,
        ) as w:
            with self.assertRaises(IndexError):
                w.write_day(self.nt, np.full((self.nlat, self.nlon), 0.5))
            with self.assertRaises(IndexError):
                w.write_day(-1, np.full((self.nlat, self.nlon), 0.5))
            # Write something so the close succeeds.
            for t in range(self.nt):
                w.write_day(t, np.full((self.nlat, self.nlon), 0.5))

    def test_rejects_wrong_shape_slab(self):
        with self.WriterCls(
            path=self.out_path,
            source_lats_deg=self.lats,
            source_lons_deg=self.lons,
            ephemeris=self.eph,
        ) as w:
            with self.assertRaises(ValueError):
                w.write_day(0, np.full((self.nlat + 1, self.nlon), 0.5))
            for t in range(self.nt):
                w.write_day(t, np.full((self.nlat, self.nlon), 0.5))

    def test_rejects_writes_after_close(self):
        w = self.WriterCls(
            path=self.out_path,
            source_lats_deg=self.lats,
            source_lons_deg=self.lons,
            ephemeris=self.eph,
        )
        for t in range(self.nt):
            w.write_day(t, np.full((self.nlat, self.nlon), 0.5))
        w.close()
        with self.assertRaises(RuntimeError):
            w.write_day(0, np.full((self.nlat, self.nlon), 0.5))

    def test_validation_rejects_bad_inputs_before_open(self):
        # Bad lats: out of range — must fail before .partial is touched.
        bad_lats = np.linspace(-100, 100, self.nlat)
        with self.assertRaises(ValueError):
            self.WriterCls(
                path=self.out_path,
                source_lats_deg=bad_lats,
                source_lons_deg=self.lons,
                ephemeris=self.eph,
            )
        # Partial file must not exist after rejection.
        self.assertFalse(self.out_path.with_name(
            self.out_path.name + ".partial").exists())

    # ---- End-to-end integration with assemble.assemble_into_writer ----

    def test_e2e_with_simple_provider(self):
        """A more practical e2e check: build a fake source provider with
        per-day-distinguishable slabs, run through ``assemble_into_writer``,
        and confirm the file holds what we sent."""
        from coordinate_transformation.assemble import assemble_into_writer

        class _PerDayValueSource:
            source_lats_deg = self.lats
            source_lons_deg = self.lons

            def slice_for_date(_self, dt):
                # Encode day-of-year fractionally so DF stays in [0, 1].
                doy = dt.timetuple().tm_yday
                v = (doy % 100) * 0.005   # 0..0.495
                return np.full((self.nlat, self.nlon), v, dtype=np.float64)

        with self.WriterCls(
            path=self.out_path,
            source_lats_deg=self.lats,
            source_lons_deg=self.lons,
            ephemeris=self.eph,
        ) as w:
            assemble_into_writer(
                _PerDayValueSource(), self.eph, w, progress_every=0,
            )

        with _nc4.Dataset(str(self.out_path), "r") as ds:
            df = ds[VAR_DF][:]
            self.assertEqual(df.shape, (self.nt, self.nlat, self.nlon))
            for t in range(self.nt):
                doy = self.eph.datetimes[t].timetuple().tm_yday
                expected = (doy % 100) * 0.005
                # Constant input → constant rotated output (rotation
                # preserves constants exactly).
                np.testing.assert_allclose(df[t], expected, atol=1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
