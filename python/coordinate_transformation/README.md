# coordinate_transformation

Post-processor that converts the Matlab output of sun-fixed dimming-factor NetCDF output from
the MATLAB PSF irradiance model into a geographic-frame NetCDF suitable for
multiplication against the CMIP6 SolarForcing TSI/SSI file, producing a CESM
input covering 2035-01-01 through 2065-12-31 (11,323 daily snapshots on the
f09 grid, 192 lat × 288 lon).

The package is entirely self-contained — a transparent, dependency-light Python
post-processor that reads the MATLAB output verbatim and emits a CF-1.8
compliant NetCDF4 file with one full global DF slab per day.

See the top-level `README.md` for the MATLAB upstream (constellation generator
and irradiance model).

---

## Scientific context

The PSF irradiance model computes the dimming factor `DF(day, lat, lon)` in a
**sun-fixed** frame: the substellar point is pinned at (lat = 0°, lon = 0°)
every day, because the model computes at noon UT with the Sun's declination
held implicit in the geometry. CESM needs the same field in a **geographic**
frame, where the Sun's declination δ(t) shifts the shadow footprint ±23.44°
with the seasons.

Each day's geographic slab is obtained by a spherical rotation of that day's
sun-fixed slab, parameterised by δ(t). The rotation is applied as a **backward
map**: for every target (φ, θ) on the output grid, we compute the source
(φ₀, θ₀) on the sun-fixed grid and bilinearly interpolate. The math is spelled
out in `rotate.py`.

δ(t) comes from JPL Horizons (geocentric declination of the Sun at 00:00 UT on
every date in the target period). Eight leap days (2036, 40, 44, 48, 52, 56,
60, 64) fall in range and are handled automatically.

---

## Pipeline

```
                                 ┌──────────────────────────┐
  Sun-fixed .nc   ──────▶ │    SunFixedSourceProvider │
  (df_standard: 365×lat×lon,     │    (slice_for_date)       │
   df_leap:     366×lat×lon)     └────────────┬──────────────┘
                                              │
                                              ▼
  JPL Horizons ephemeris ─────▶  ┌──────────────────────────┐
  (date strings, δ[deg])         │  assemble_into_writer    │
                                 │  (for each day t:        │
                                 │     src  = source[t]     │
                                 │     dst  = rotate(src,δ) │──┐
                                 │     write_day(t, dst))   │  │
                                 └──────────────────────────┘  │
                                              │                │
                                              ▼                │
                                 ┌──────────────────────────┐  │
                                 │    CesmNcWriter          │◀─┘
                                 │    (CF-1.8, chunked,     │
                                 │     compressed, atomic)  │
                                 └────────────┬─────────────┘
                                              ▼
                               df_2035_2065.nc   (geographic, CF-compliant)
```

---

## Modules

**`horizons_parser.py`** — Reads a JPL Horizons "text" ephemeris file
(`$$SOE` / `$$EOE` delimited), parses the date and declination on each row,
converts Horizons' "±DD MM SS.s" sign convention into signed decimal degrees,
and returns a `HorizonsEphemeris` dataclass (`date_strings`, `datetimes`,
`subsolar_lat_deg`). Parsing is strict: bad rows raise with row numbers so the
user can repair the file.

**`rotate.py`** — Pure numpy implementation of Illeana's spherical rotation
routine. `rotate_shade_field(src, src_lats, src_lons, out_lats, out_lons, δ)`
returns the rotated 2-D slab. Helpers: `backward_map` (the analytic pullback
from output to source coordinates) and `_bilinear_interp_with_lon_wrap` (a
bilinear interpolator that pads the source grid with one wrap-around longitude
column so queries near 360°/0° are periodic by construction). A
`RotationDiagnostics` return path is available for audit runs
(`return_diagnostics=True`).

**`assemble.py`** — The orchestrator. It defines two protocols — a
`SourceProviderProtocol` that returns a daily slab given a `datetime`, and a
`WriterProtocol` that accepts hyperslab writes — and a pure-Python
`assemble_into_writer` function that iterates over the ephemeris dates,
selects the right sun-fixed slab (standard vs. leap calendar), rotates it,
and streams the result to the writer. The top-level `assemble_from_files`
glues the real I/O together. A `--show-progress` path pulls in `tqdm` lazily
and falls back to periodic log lines if `tqdm` is not installed.

**`write_cesm_nc.py`** — The NetCDF writer. `CesmNcWriter` is a context
manager: on entry it opens `<output>.partial`, writes the full header
(dimensions, coordinate variables with CF attributes, pre-allocated `DF`
chunked at `(1, nlat, nlon)` with zlib+shuffle at level 4), and pre-allocates
the `date(time, 10)` character array with ISO 8601 strings from the ephemeris.
Each `write_day(t, slab)` is a single hyperslab write. On clean close, the
partial file is atomically renamed to the final path; on exception inside the
`with` block, the partial file is removed. The final file therefore either
exists fully or not at all — no half-written 5 GB NetCDF after a crash.

---

## Installation

```
pip install -r python/coordinate_transformation/requirements.txt
```

Runtime requirements are minimal: `numpy` and `netCDF4`. A progress bar is
available if `tqdm` is installed (optional — the code logs periodic progress
without it).

Tests use Python's stdlib `unittest` only; no separate test-runner install is
required.

---

## Usage

### Programmatic (recommended for reproducible runs)

```python
from coordinate_transformation.assemble import assemble_from_files

assemble_from_files(
    source_nc_path   = "path/to/source.nc",
    horizons_txt_path = "path/to/horizons.txt",
    output_nc_path   = "df_2035_2065.nc",
    show_progress    = True,          # tqdm bar if installed, logs if not
)
```

### Script invocation

Edit the paths in `coordinate_transformation/transform.py` to point to your
source NetCDF, Horizons ephemeris, and desired output path, then run from the
`python/` directory:

```
cd path/to/python
python -m coordinate_transformation.transform
```

The `-m` flag is required (running the file directly with `python
coordinate_transformation/transform.py` will fail with a `ModuleNotFoundError`
because Python won't add the parent directory to `sys.path`).

---

## Input contract

### Source NetCDF 

Must contain two DF fields covering a standard and a leap year:

| variable              | shape               | meaning                                   |
|-----------------------|---------------------|-------------------------------------------|
| `dimming_factor_standard` (or `df`)  | `(365, nlat, nlon)` | one sun-fixed slab per DOY in a 365-day year |
| `dimming_factor_leap`     (or `dfl`) | `(366, nlat, nlon)` | one sun-fixed slab per DOY in a 366-day year |

Lat/lon coordinate variables must be 1-D, strictly increasing, in degrees,
covering [-90, 90] and [0, 360) respectively. If the source file's lat/lon
arrays are stubbed (e.g., integer template indices), the CESM f09 grid is
used as a fallback and a warning is logged.

The time axis may be first (`(t, lat, lon)`) or last (`(lat, lon, t)`); the
reader auto-reorients using a heuristic (matching 365/366 against the axis
lengths). Ambiguous cases (e.g., `365 × 365 × 7`) are rejected with a clear
error.

### Horizons ephemeris

Plain-text JPL Horizons output with `$$SOE` / `$$EOE` markers. The parser
reads date (col 1, `YYYY-Mmm-DD` format) and declination (col 3 or 4, in
`±DD MM SS.s` sexagesimal). One row per day at 00:00 UT. Must span the full
target period inclusive of both endpoints.

---

## Output schema

```
dimensions:
    time       = 11323               // days, 2035-01-01 through 2065-12-31
    latitude   = 192                 // f09
    longitude  = 288                 // f09
    string_len = 10                  // len("YYYY-MM-DD")

variables:
    double   time      (time)                     "days since 2035-01-01 12:00:00" / gregorian
    double   latitude  (latitude)                 degrees_north
    double   longitude (longitude)                degrees_east
    char     date      (time, string_len)         ISO 8601 calendar date
    double   DF        (time, latitude, longitude)
                                                   long_name = "dimming factor"
                                                   units     = "1"
                                                   valid_range = 0.0, 1.0
                                                   _FillValue  = 9.96921e36
                                                   chunking = (1, 192, 288), zlib+shuffle, level 4
```

Global attributes: `Conventions = "CF-1.8"`, plus `title`, `institution`
(Planetary Sunshade Institute), `source` (names the rotation module),
`references` (GitHub repo URL), `history` (UTC creation timestamp + input
paths + period), and a `comment` explaining that DF is intended to be
multiplied externally against the CMIP6 SolarForcing TSI/SSI.

Variable names follow Illeana's spec (`latitude` / `longitude`), not the
shorter `lat`/`lon`, because this file feeds her external multiplier.

---

## Transparency features

The project mandate is "lots of visibility into what is happening … transparent
to a wide scientific audience." Concretely:

- **Day-0 audit log**: on the first day processed, `assemble_into_writer` logs
  the date, declination, source-slab statistics (min/mean/max), and output-slab
  statistics so an operator can sanity-check the pipeline from a single log
  line before committing to an 11-thousand-day run.
- **Plausibility check**: every rotated slab is checked against DF ∈ [0, 1]
  (with 1e-6 slack for floating-point). Out-of-range values halt the run with
  the date printed.
- **Progress reporting**: `tqdm` bar if installed; otherwise periodic log lines
  at a configurable interval.
- **Self-describing output**: `ncdump -h df_2035_2065.nc` shows the GitHub URL,
  the source files used, the creation timestamp, and the period covered. A
  scientist with only the `.nc` file can retrace its origins.
- **Atomic file writes**: the final file either exists completely or not at
  all. No half-written 5 GB NetCDFs after a crash.
- **Pure-numpy rotation**: no scipy interpolator black box — the bilinear
  interpolation and backward map are readable, ~40 lines of numpy, and carry
  the analytic formulas in comments.
- **Dataclass return types and Protocol-based seams**: the orchestrator's
  collaborators are all explicit, typed, and mockable — no hidden I/O in the
  pure code paths, and every module has dependency-free unit tests.

---

## Running the tests

```
cd (to source directory))
python3 -m unittest discover -s python/coordinate_transformation/tests -v
```

Full suite: **117 tests** across four modules. Sixteen of those are writer
integration tests that exercise real NetCDF I/O; they skip cleanly in
environments without `netCDF4` installed. In a fully-provisioned environment
all 117 run green.

Each module has its own test file:

| module            | test file                                | tests  | notes                                        |
|-------------------|------------------------------------------|--------|----------------------------------------------|
| `horizons_parser` | `tests/test_horizons_parser.py`          | ~30    | incl. real-file regression (11,323 entries)  |
| `rotate`          | `tests/test_rotate.py`                   | 23     | identity, solstice geometry, area-weighted invariance, f09 grid |
| `assemble`        | `tests/test_assemble.py`                 | 32     | mocks writers via Protocol; real ephemeris integration          |
| `write_cesm_nc`   | `tests/test_write_cesm_nc.py`            | 31     | 15 pure helpers + 16 NetCDF integration (auto-skip w/o netCDF4) |

---

## Design notes

**Pure numpy over scipy.** The rotation is implemented in ~40 lines of numpy
(backward map + bilinear interp with wraparound) rather than via
`scipy.RegularGridInterpolator`, because the mandate is transparency: every
scientist reading the code should be able to trace every arithmetic step to
the derivation in `rotate.py`'s docstring.

**Protocol-based dependency injection.** `assemble.py` defines
`SourceProviderProtocol` and `WriterProtocol` so the orchestrator can be tested
with mock writers and mock sources — no netCDF4, no filesystem — while
production I/O is done through thin glue classes (`SunFixedSourceProvider`,
`CesmNcWriter`) that satisfy those protocols.

**Area-weighted means, not naive means.** A constant input DF must round-trip
to a constant output under rotation. Tests verify this using `cos(lat)`-area-
weighted means, because the f09 grid's polar cells are narrower in area than
equatorial cells — the naive `np.mean` over-weights the poles and the test
would falsely fail.

**Longitude wrap.** The source grid is padded with one wrap-around column
(copy of lon=0 appended at lon=360) before bilinear interpolation, so queries
near the seam are periodic by construction. This avoids a separate branch in
the interpolator.

**Atomic write via `.partial` rename.** On a 5 GB output file, an unhandled
exception partway through would otherwise leave a corrupt file that looks
superficially like a valid NetCDF. Writing to `<path>.partial` and renaming
on clean close (and unlinking on exception) gives us "all or nothing"
persistence.

**Lazy imports for optional deps.** `netCDF4` is imported inside
`CesmNcWriter.__init__`; `tqdm` is imported inside the progress code path.
This keeps pure-logic tests runnable in a minimal sandbox without scientific
stacks.

---

## Known limitations

- **CESM rotation logic still pending.** CESM's `solar_shade.F90` may re-rotate
  the field when applying it in a climate run; the PSI / CESM science team
  should confirm that the geographic frame we emit here is what
  `solar_shade.F90` expects. If an additional azimuthal rotation is needed to
  match CESM's conventions, it is cleanest to apply it here in Python (one
  well-tested rotation) rather than in F90.
- **Polar interpolation.** Near the poles the bilinear interpolation on a
  lat/lon grid is less accurate than spherical triangle interpolation, but
  the source model itself is grid-aligned and the dimming factor is smooth
  there, so we accept the bilinear error in exchange for transparency.
- **f09 only, tested.** The pipeline is written to accept any regular lat/lon
  grid for source and output (and validates non-monotonic / out-of-range
  inputs), but only the f09 grid (192 × 288) has been exercised end-to-end.
- **Single-year source repeated across 31 years.** Upstream MATLAB
  model currently produces a single representative year (standard + leap);
  this post-processor replays that year against each calendar year in
  [2035, 2065]. Interannual variation in the constellation geometry is a
  future-work item (see the top-level README, "Grow constellation over time").
