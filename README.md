# Planetary Sunshade Irradiance Model

**Planetary Sunshade Institute — v26.04.07**

This code is written to support planetary sunshade researchers in generating a constellation of heliogyro sunshades, calculating the shade on the Earth for that array, and outputting a NetCDF (.nc) file which can be read by climate models. Given a described array of sunshade elements near the Sun-Earth L1 point, the model computes the irradiance factor across a latitude-longitude grid of the Earth's surface and exports the result compatible with the NCAR Community Earth System Model (CESM).

---

## Mathematical Foundation

The model is built on the mathematical framework described in:

**Raymond, M. (2026). *Irradiance + Sunshade Modeling*, v26.03.18. Planetary Sunshade Foundation.**

This document, included in the repository as `v26_03_18_pfl_sunshade_model.pdf`, defines the complete mathematical basis for the model. It establishes a general surface partition framework, derives the irradiance and shaded irradiance approximations, and describes the limb darkening treatment. All subsequent development — including the bug fixes and CESM interface work documented below — is built on this foundation. Researchers using or extending this code should read the PDF first.

---

## Overview

The model first requires a constellation to be defined and generated. The Constellation Generator places heliogyro sunshades in an array according to user-specified parameters, producing a `.mat` file that the irradiance model reads as input.

The irradiance model then propagates solar flux from the Sun to the Earth through that intermediate layer of sunshade elements positioned near the Sun-Earth L1 Lagrange point. For each partitioned subsurface of the Earth, it can compute:

- **Irradiance**: total solar flux incident on that subsurface
- **Shaded irradiance**: flux blocked by the shade constellation
- **Shading factor**: ratio of shaded to total irradiance
- **Irradiance factor**: 1 minus the shading factor — the quantity written to the CESM NetCDF file

The model accounts for solar limb darkening, the 1/r² scaling of flux, the geometric projection of shade elements onto the solar disk, and the time-varying Sun-Earth geometry via SPICE-derived kinematics data.

---

## Repository Structure

```
matlab/sunshade models/psf/current version/   — MATLAB source code
excel/psf model/                              — Excel and .mat input files; constellation outputs
numerical control/                            — .nc reference and template files
Reports/                                      — Investigations into functionality
v26_03_18_pfl_sunshade_model.pdf              — Mathematical documentation
config_paths.m                                — Portable path configuration (edit once per machine)
```

---

## Getting Started

### 1. Path setup

Paths are managed centrally through `config_paths.m` at the repository root. On a fresh clone, no configuration is needed — `config_paths()` detects the project root automatically from its own file location. If you have multiple clones or a non-standard setup, set the environment variable `SUNSHADE_PROJECT_ROOT` to your repo path (e.g. in your MATLAB `startup.m`).

Add the `current version` folder and all its subfolders to the MATLAB path.

### 2. Generate a constellation

Before running the irradiance model, you need a constellation file. Open and edit `constellation generator/input_Constellation_V2_Parameters.m` to define your array. Key parameters include cluster geometry, placement method (`lattice` or `random`), and any time-varying motion function. Then run:

```matlab
generate_constellation_v2
```

This produces a `.mat` file in `excel/psf model/`, named with the architecture label, craft count, and date. It also produces a diagnostic image.

### 3. Configure the irradiance model

Open `input/input_Model_Parameters___0___Sr.m` and set:

- **Analysis type**: typically `'irradiance factor'`
- **Configuration**: `'manual'` for a quick preview, or `'preconfigured: low resolution'` to generate a full CESM-compatible NetCDF
- **Time parameters**: start date, frequency (days between outputs), and number of periods

Open 'file locations/location_Heliogyro_Kinematics_Data___E___Sr.m'
- **Constellation file**: point the excel_file_name to the `.mat` file generated in step 2

### 4. Preview in manual mode

With `configuration = 'manual'`, run `analysis_General___0___X` from the MATLAB command window. This outputs a small irradiance matrix to the command window and opens a diagnostic chart — a latitude × month map of the irradiance factor — so you can verify the shade pattern is physically reasonable before running a full year.

### 5. Generate the NetCDF output

Switch to `configuration = 'preconfigured: low resolution'` and re-run `analysis_General___0___X`. This iterates over the full year and a leap year at your specified frequency and writes the 192-latitude × 288-longitude × time irradiance factor grid to a NetCDF file in `numerical control/exports/`. It also saves the diagnostic chart as a .png alongside the .nc file.

---

## Constellation Generator V2

### Overview

The V2 generator creates time-varying constellations of heliogyro sunshades and exports them as `.mat` files for use by the irradiance model. Constellations are made up of one or more named **clusters**, each occupying an elliptical footprint in the Y-Z plane at a defined position in L1-centred synodic coordinates. Each cluster can move over the course of the year via a user-supplied motion function.

### Input file

All parameters are set in `input_Constellation_V2_Parameters.m`. The file is structured and commented for direct editing; no other files need to be changed for typical runs.

### Placement methods

**Lattice** (recommended): Hexagonal close-packing with FCC 3D offset across planes. Craft count is determined by geometry — the lattice fills the ellipse efficiently with guaranteed zero line-of-sight overlap. No minimum-buffer enforcement is needed.

**Random**: Craft positions are sampled randomly within the ellipse boundary, with a minimum edge-to-edge buffer enforced between sails. Craft count is set explicitly by `n_craft`.

### Coordinate convention

All positions are in km, in the L1-centred synodic frame. X positive = sunward (toward the Sun), Y = tangential to Earth's orbital motion, Z = ecliptic north. A positive X value places a shade element sunward of L1, which is necessary to offset solar radiation pressure on a reflective sail. The optimal shading equilibrium point is approximately 838,400 km sunward of L1 (McInnes, 2002).

### Heliogyro position limits and first-order validation

The heliogyro position file defines the allowable region in which a heliogyro can fly, based on rigid blade dynamics and a fixed normal pointing toward the Sun. This region has been found to be more than ample for Earth-relevant shading geometries. Craft positions are checked against this boundary as a first-order validation step during constellation generation.

### Normal vectors

The normal vector columns (NX, NY, NZ) in the kinematics file support future analysis in which individual heliogyros may be tilted relative to the Sun-Earth axis — for example, to steer shade toward specific latitudes or to implement out-of-plane orbital configurations of the type studied by Sánchez & McInnes (2015). This capability is implemented in the model but not yet exercised in current array designs.

### Output

The generator writes a `.mat` file to `excel/psf model/`, named with the architecture label, craft count, and date. A preview image of the constellation geometry (face-on Y-Z projection) is saved alongside it. A diagnostic chart opens at run time showing the constellation layout.

### Spectral data

`star spectral data.xlsx` defines the solar spectrum and limb darkening coefficients used in the irradiance calculation. The default values (C1 = 0.61, C2 = 0.39) implement the Eddington linear limb darkening approximation (Tripathi et al., 2020). Multiple spectral bands can be defined; the model produces a separate irradiance result per band.

---

## File Development History

### v26.03.18 → v26.03.19: Initial validation and bug fixes

During validation against Sánchez & McInnes (2015), three bugs were identified and corrected. See `Reports/Sanchez-McInes static array validation analysis.docx` for the full validation methodology and results.

**Bug 1 — Incorrect L1 distance calculation** (`distance_PL1___Sr___Sr.m`): `fzero` was converging to a spurious root of the L1 equilibrium equation, placing shade elements ~227,000 km too close to Earth. Fixed by replacing with the standard first-order approximation: `r_L1 = r_SP × (M_earth / 3M_sun)^(1/3)`, accurate to < 0.01%.

**Bug 2 — Dot product computed as element-wise multiply** (`analysis_Shaded_Irradiance___Sr___nM.m`): The cosine projection of sail area was computed with `.*` instead of `sum(..., 1)`, producing a 3×N matrix rather than the correct 1×N scalar per sail. This caused approximately half the expected shading to be lost.

**Bug 3 — Hardcoded Excel range limiting sail count to 11** (`location_Heliogyro_Kinematics_Data___E___Sr.m`): Ranges extended to row 20,000; blank rows are filtered downstream by existing NaN-removal logic.

Following these fixes, the model returned a center-point shading factor of **1.96%** against the Sánchez & McInnes published value of ~1.9% (3.2% error). Penumbra gradient, full-disk coverage, and latitudinal symmetry were all correctly reproduced.

---

### v26.03.19 → v26.03.25: Limb darkening investigation

A detailed investigation established that limb darkening is correctly implemented and physically consistent. Key findings: the Eddington linear approximation (I/I₀ = 0.61 + 0.39·cosθ) is appropriate for total solar intensity; the crossover point at which limb darkening makes a shade more vs. less effective than a uniform disk falls at r/R ≈ 0.83; and concentrating the constellation near the optical axis (as most practical designs do) produces slightly higher shading efficiency than a uniform-disk model would predict. See `Reports/limb_darkening_walkthrough.md` and `Reports/limb_darkening_implications.html`.

---

### v26.03.25 → v26.04.03: Shading gradient bug fixes

Analysis of constellation runs revealed a spurious ~38–40% dropoff in shading factor from the equator to the poles — far exceeding the ~7–9% expected from limb darkening alone. Two bugs were identified and fixed. See `Reports/shading_gradient_bug_fix_notes.md` for the full analysis.

**Bug 4 — cos²(θ_S) double-counting** (`analysis_Shaded_Irradiance___Sr___nM.m`): The solar disk projection factor cos(θ_S) was applied twice — once when projecting shade areas and again in the limb darkening term — producing a spurious cos²(θ_S) dependence responsible for ~31 percentage points of the false gradient.

**Bug 5 — Asymmetric limb darkening** (`analysis_Irradiance___Sr___nM.m`): Limb darkening had been inadvertently disabled in the irradiance (denominator) function during earlier testing but remained active in the shaded irradiance (numerator). This caused a ~3× underestimate of absolute shading factor at all latitudes, without affecting the gradient. Both functions now apply limb darkening identically.

After both fixes, the equator-to-pole gradient is 7.32%, consistent with the expected limb darkening effect.

---

### v26.04.03 → v26.04.07: CESM interface investigation and axial tilt fix

A systematic audit of the complete PSF–CESM interface was conducted, accounting for every physical source of solar irradiance variation and establishing clear ownership (PSF vs. CESM-internal) for each. The audit is documented in full in `Reports/CESM and MATLAB Audit.xlsx`. Supporting technical notes covering orbital geometry, solar source properties, coordinate systems, and the CESM `solar_shade` pipeline are in `Reports/`.

The primary finding was a coordinate mapping error in the PSF output:

**Bug 6 — Axial tilt coordinate mapping** (`analysis_General___0___X.m`): The PSF code computes the irradiance factor in a heliocentric, Earth-centred frame where φ = 0 corresponds to the sub-stellar point, not geographic latitude = 0°. Previously, the 192-row result was written directly to the NetCDF file as geographic latitude, pinning the shadow footprint at the equator year-round regardless of season — an error of up to ±23.44°.

The fix implements solar declination δ(day) via the Spencer (1971) Fourier series (new file: `spencer_declination.m`, accurate to ±0.035°). A 27-cell buffer is appended to both poles of the 192-row result before each daily write, and the 192-row output window is extracted at an offset of `round(δ / 0.9424°)` rows. This correctly shifts the shadow footprint to the geographic latitude of the sub-stellar point for each day of year, with an integer rounding error ≤ ±0.47° (half a cell).

Verification (`test_declination_shift_verification.m`): the shadow trough (minimum irradiance factor) tracks the expected sub-stellar latitude on all five test dates (January 1, March 21, June 21, September 22, December 21) to within the rounding margin. See `Reports/cesm_axial_tilt_report.docx` for the full analysis.


---

## Future Work and Known Issues

### Active development

**Longitude averaging and grid representation** — The current implementation takes a single latitude strip and replicates it uniformly across all 288 CESM longitude grid points. This is a  correct approximation for a constellation symmetric around the Earth-Sun axis, but the slight decrease in shade across the radius of the earth will decrease the average shade on the middle strip slightly. Plus, future constellations may have longitudinally-referenced shading. 

**L1 distance reconciliation** — Reconcile the current first-order L1 approximation with Matt's updated orbital mechanics work.

**Blade Angle θ** — Implement the ability for shades to vary their overall blade angle θ and adjust their distance from the Sun accordingly. A shade tilted away from normal moves to a new equilibrium point further from the Sun, reducing its effective shade area but allowing preferential shading of one hemisphere. The intended design: maximum shading intensity during northern hemisphere summer, reduced but still significant shading during southern hemisphere summer, with the constellation "in motion" between the two. Sub-task: establish the correlation between distance from the Sun and average blade angle θ using reference data.

**Genetic algorithm constellation optimization** — Explore genetic algorithms as an approach to optimizing constellation geometry for target shade outcomes.

**Grow constellation over time** - Currently we repeat the same year's worth of data in the model. How do we model a ramp up? 

### Lower priority cleanup

1. ~~Generalize filepath configuration~~ *(completed — see `config_paths.m`)*
2. Project L1 position / SPICE data forward through all future dates.
3. Split `material_irradiance_absorption` into two functions: one for total geometric area, one for material opacity (1.0 for heliogyro).
4. Incorporate the effect of aphelion/perihelion motion on L1 equilibrium dynamics throughout the year.

---

## Other Notes

`L1_Stability_Region_Data_Visualization.m` reads the heliogyro kinematics file and produces diagnostic plots of the constellation geometry: a 3D point cloud, orthogonal projections (XY, XZ, YZ), position histograms, and an outer envelope plot of the array face-on. Useful for spot-checking constellation geometry outside of a full model run. Set the `filename` path in the User Settings section at the top, then run directly in MATLAB.

---

## References

- Raymond, M. (2026). *Irradiance + Sunshade Modeling*, v26.03.18. Planetary Sunshade Foundation.
- Sánchez, J.-P. and McInnes, C.R. (2015). Optimal Sunshade Configurations for Space-Based Geoengineering near the Sun-Earth L1 Point. *PLOS ONE*, 10(8): e0136648. https://doi.org/10.1371/journal.pone.0136648
- McInnes, C.R. (2002). Minimum Mass Solar Shield for terrestrial climate control
- Spencer, J.W. (1971). Fourier series representation of the position of the sun. *Search*, 2(5), 172.
- Tripathi, D. et al. (2020). Study of Limb Darkening Effect and Rotation Period of Sun by using Solar Telescope. *Journal of Scientific Research*, 64(1).

---

*Planetary Sunshade Institute — planetarysunshadeinstitute on GitHub*
