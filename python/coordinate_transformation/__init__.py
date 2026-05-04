"""
coordinate_transformation
=========================

Post-processor that converts Matt's sun-fixed dimming-factor NetCDF output into
a geographic-frame NetCDF suitable for multiplication against the CMIP6
SolarForcing file.

Scientific context
------------------
The PSF irradiance model (MATLAB) produces a dimming factor field DF(day, lat,
lon) in a sun-fixed frame: the substellar point sits at (lat=0, lon=0) on each
day's grid, because the model computes at noon UTC with the Sun's declination
treated implicitly in the shade geometry. For CESM ingestion we need the same
field expressed in geographic (geographic lat, geographic lon) coordinates,
where the Sun's declination δ(t) shifts the shadow footprint ±23.44° with the
seasons.

This package reads:
  - Matt's existing 2-year .nc file (variables `dimming_factor_standard` and
    `dimming_factor_leap`, or equivalently `df` / `dfl`; they're aliased in the
    template to the same underlying datasets).
  - The JPL Horizons ephemeris file giving the Sun's geocentric declination
    δ(t) at 00:00 UT for every day in the target period (2035-01-01 through
    2065-12-31).

...and writes a single combined DF(time, latitude, longitude) field where each
day's sun-fixed slab has been rotated by that day's δ. The output carries
CF-compliant time and calendar attributes so downstream preprocessing tools can
multiply it element-wise against CMIP6's SSI(time, wavelength) via numpy/xarray
broadcasting.

The rotation algorithm follows Illeana's Python specification (see the RTF
instruction file in `Coordinate Transformation/`). See `rotate.py` for the
full derivation.

Modules
-------
horizons_parser : read the JPL Horizons text file, return (date, δ) per day.
rotate          : rotate a 2D sun-fixed shade field into geographic coords.
assemble        : orchestrator — glues parser, rotation, and writer together.
write_cesm_nc   : NetCDF writer (dims, vars, CF attributes).
"""
