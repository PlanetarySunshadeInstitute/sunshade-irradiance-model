# PSF Irradiance and Sunshade Model

**Planetary Sunshade Foundation — v26.03.19**

This code is written to support planetary sunshade researchers in taking the input of a constellation of sunshades and outputting a file which can be read by climate models. Given a described array of sunshade elements near the Sun-Earth L1 point, the model computes the irradiance factor across a latitude-longitude grid of the Earth's surface and exports the result as a NetCDF file compatible with the NCAR Community Earth System Model (CESM).

---

## Mathematical Foundation

The model is built on the mathematical framework described in:

**Raymond, M. (2026). *Irradiance + Sunshade Modeling*, v26.03.18. Planetary Sunshade Foundation.**

This document, included in the repository as `v26_03_18_pfl_sunshade_model.pdf`, defines the complete mathematical basis for the model. It establishes a general surface partition framework, derives the irradiance and shaded irradiance approximations, and describes the limb darkening treatment. All subsequent development — including the bug fixes documented below — is built on this foundation. Researchers using or extending this code should read the PDF first.

---

## Overview

The model propagates solar flux from the Sun to the Earth through an intermediate layer of sunshade elements positioned near the Sun-Earth L1 Lagrange point. For each partitioned subsurface of the Earth, it computes:

- **Irradiance**: total solar flux incident on that subsurface
- **Shaded irradiance**: flux blocked by the shade constellation
- **Shading factor**: ratio of shaded to total irradiance
- **Irradiance factor**: 1 minus the shading factor — the quantity written to the CESM NetCDF file

The model accounts for solar limb darkening, the 1/r² scaling of flux, the geometric projection of shade elements onto the solar disk, and the time-varying Sun-Earth geometry via SPICE-derived kinematics data.

---

## Repository Structure

```
sunshade models/psf/current version/   — MATLAB source code
excel/psf model/                       — Excel input files
numerical control		       — .nc reference and template files
v26_03_18_pfl_sunshade_model.pdf       — Mathematical documentation
```

---

## Getting Started

See Section 2 of the documentation PDF for full setup instructions. In brief:

1. Add the `current version` folder and its subfolders to the MATLAB path
2. Set the paths to your `excel` and `numerical control` folders in `location_Excel_Folder___0___St.m` and `location_NC_Folder___0___St.m`
3. Configure analysis parameters in `input/input_Model_Parameters___0___Sr.m`
4. Run `analysis_General___0___X` from the MATLAB command window

A successful test run outputs an 11×11 matrix of irradiance values to the command window.

---

## Excel Input Files

The model reads two categories of Excel input:

### Heliogyro kinematics data (`heliogyro kinematics data.xlsx`)

This file defines the position and orientation of each sunshade element in the constellation in relationship to Sun Earth Lagrange 1 (L1). Each row represents one heliogyro unit with six columns:

| Column | Description | Status |
|--------|-------------|--------|
| PX | X position in km, relative to L1 (along Sun-Earth axis) | Active |
| PY | Y position in km, relative to L1 (tangential to planet motion) | Active |
| PZ | Z position in km, relative to L1 (normal to orbital plane) | Active |
| NX | X component of surface normal vector | Reserved for future use |
| NY | Y component of surface normal vector | Reserved for future use |
| NZ | Z component of surface normal vector | Reserved for future use |

A positive value for the X position describes a sunshade element sunward of L1, which is necessary to offset the solar radiation pressure on a sail. The optimal shade-to-mass equilibrium point for sunshades is considered to be 2.36 million km from Earth.

The normal vector columns (NX, NY, NZ) are included to support future analysis in which individual heliogyros may be tilted relative to the Sun-Earth axis — for example, to steer shade toward specific latitudes or to implement out-of-plane orbital configurations of the type studied by Sánchez & McInnes (2015). This capability is implemented in the model but not yet exercised in current array designs.

### Star spectral data (`star spectral data.xlsx`)

This file defines the solar spectrum and limb darkening coefficients used in the irradiance calculation. The default values are cited from Tripathi et al. (2020) for total intensity. Multiple spectral bands can be defined; the model outputs a separate irradiance result for each band.

---

## Diagnostic Tools

Two visualization scripts are included to support analysis and verification:

### `plot_shading_map.m`

Reads a completed NetCDF output file and renders a world map showing the irradiance factor at a specified time step. The map includes coastline overlay, a color scale from full sun to fully shaded, and summary statistics (mean irradiance factor, minimum irradiance factor, fraction of Earth shaded). Useful for verifying that the shade pattern is physically reasonable before submitting to CESM.

To use: set the `nc_file` path and `time_index` in the User Settings section at the top of the script, then run directly in MATLAB.

### `L1_Stability_Region_Data_Visualization.m`

Reads the heliogyro kinematics Excel file and produces a set of diagnostic plots showing the spatial distribution of the sail constellation: a 3D point cloud, orthogonal projections (XY, XZ, and YZ), histograms of the position distributions, and an outer envelope boundary plot of the array face-on. This allows researchers to visually inspect the constellation geometry — including array footprint, radial spread, and depth along the Sun-Earth axis — before running the full irradiance calculation.

To use: set the `filename` path in the User Settings section, then run directly in MATLAB.

---

## Bug Fixes: v26.03.18 → v26.03.19

During validation of the model against Sánchez & McInnes (2015), three bugs were identified and corrected. These are documented here in full, as they are non-trivial and affect the physical correctness of the results.

### Bug 1: Incorrect L1 distance calculation

**File**: `operations/distance_PL1___Sr___Sr.m`

**Symptom**: The model placed shade elements approximately 227,000 km too close to Earth, causing the shade's projected area on the Sun to be dramatically underestimated. Most planet surface points saw zero sail intersections with the Sun.

**Root cause**: The L1 equilibrium equation was solved using `fzero` with an initial guess of `0.01 × R`. The equation has two roots, and `fzero` was consistently converging to a spurious root at ~1,294,208 km rather than the physical L1 distance of ~1,521,504 km. Evaluation of the equation at both roots confirmed that neither was exactly zero — the simplified three-body equilibrium equation as implemented does not accurately reproduce the full quintic equation governing the collinear Lagrange point.

**Fix**: Replaced the `fzero` call with the standard first-order approximation for the L1 distance:

```matlab
% Before (buggy):
kinematics_SP_data.distance.PL1 = fzero(L1_equilibrium_equation, 0.01 * R);

% After (fixed):
kinematics_SP_data.distance.PL1 = kinematics_SP_data.distance.SP * (masses.planet / (3 * masses.star))^(1/3);
```

This approximation gives 1,521,600 km, matching the SPICE-derived expected value of 1,521,504 km to within 100 km (< 0.01% error).

---

### Bug 2: Projection computed as element-wise multiply instead of dot product

**File**: `analyses/analysis_Shaded_Irradiance___Sr___nM.m`

**Symptom**: Even after fixing bug 1, the global average shading was approximately half the expected value.

**Root cause**: The cosine projection of each sail's area onto the planet-to-star direction was computed using MATLAB's element-wise multiply operator (`.*`) on two 3×N matrices, producing a 3×N result rather than the intended 1×N scalar per sail. A dot product requires multiplying corresponding components and **summing across the three spatial dimensions**, collapsing the result to a single scalar per sail representing the cosine of the angle between the sail's normal and the viewing direction.

**Fix**:

```matlab
% Before (buggy) — produces a 3×N matrix:
projected_shade_areas = area_effective * vectors.unit.p.S .* unit_normal_vectors_S;

% After (fixed) — produces a 1×N scalar per sail:
projected_shade_areas = area_effective * sum(vectors.unit.p.S .* unit_normal_vectors_S, 1);
```

---

### Bug 3: Hardcoded Excel range limiting sail count to 11

**File**: `file locations/location_Heliogyro_Kinematics_Data___E___Sr.m`

**Symptom**: Only 11 sails were loaded regardless of how many were present in the Excel file.

**Root cause**: The Excel ranges for reading heliogyro positions and normal vectors were hardcoded to rows 11–21, loading exactly 11 rows.

**Fix**: Extended the ranges to row 10,000, accommodating arrays of up to ~10,000 heliogyros. Blank rows beyond the last sail are automatically filtered downstream by existing NaN-removal logic.

```matlab
% Before (buggy):
excel_file.ranges.position_vectors = 'A11:C21';
excel_file.ranges.normal_vectors   = 'D11:F21';

% After (fixed):
excel_file.ranges.position_vectors = 'A11:C10000';
excel_file.ranges.normal_vectors   = 'D11:F10000';
```

---

## Validation

Following the bug fixes, the model was validated against the baseline case from:

**Sánchez, J.-P. and McInnes, C.R. (2015). Optimal Sunshade Configurations for Space-Based Geoengineering near the Sun-Earth L1 Point. *PLOS ONE*, 10(8): e0136648.**

A single solid opaque disk of 1,434 km radius placed 2.44 × 10⁶ km from Earth was modeled. The PSF model returned a center-point shading factor of **1.96%**, compared to the published value of approximately **1.9%** — an error of 3.2%. The penumbra gradient, full-disk coverage, and latitudinal symmetry were all correctly reproduced. This confirms that the model's shading geometry, solar limb darkening implementation, and flux propagation are physically correct.

A full account of the validation methodology and results is provided in the companion document `validation_section.docx`.

---

## References

- Raymond, M. (2026). *Irradiance + Sunshade Modeling*, v26.03.18. Planetary Sunshade Foundation.
- Sánchez, J.-P. and McInnes, C.R. (2015). Optimal Sunshade Configurations for Space-Based Geoengineering near the Sun-Earth L1 Point. *PLOS ONE*, 10(8): e0136648. https://doi.org/10.1371/journal.pone.0136648
- Tripathi, D. et al. (2020). Study of Limb Darkening Effect and Rotation Period of Sun by using Solar Telescope. *Journal of Scientific Research*, 64(1).

---

*Planetary Sunshade Foundation — planetarysunshadeinstitute on GitHub*
