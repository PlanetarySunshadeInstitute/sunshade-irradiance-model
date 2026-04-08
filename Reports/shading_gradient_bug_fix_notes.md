# Shading Gradient Bug Fix Notes
**Date:** April 3, 2026
**Files changed:** `analysis_Shaded_Irradiance___Sr___nM.m`, `analysis_Irradiance___Sr___nM.m`

---

## Observed Problem

Runs of `analysis_General` with various constellations, using both `disc` and `sphere` planet partition types, showed an ~38–40% dropoff in `shaded irradiance`, `shading factor`, and `irradiance factor` when comparing the equator to the poles of the output grid. The expected dropoff is ~7–9%, arising from limb darkening (the sun's disc is dimmer near the edges, so a shade shifted toward the limb blocks less light).

---

## Two Bugs Found

### Bug 1 — cos²(θ_S) double-counting in `analysis_Shaded_Irradiance`
**Responsible for ~31 percentage points of the false gradient.**

**The physics.** A shade element near L1 casts a shadow cross-section *A* onto the sun's surface. At a point on the sun's disc at angle θ_S from the disc normal, the shadow footprint on the actual (tilted) surface is *A / cos(θ_S)*. The power emitted by that footprint toward the observer follows Lambert's cosine law: *P ∝ B · cos(θ_S) · A_surface*, where *B* is the surface brightness. These two factors cancel:

```
P ∝ B · cos(θ_S) · (A / cos(θ_S)) = B · A
```

The blocked irradiance is **independent of where on the disc the shadow falls**. The sun's disc appears uniformly bright (before limb darkening), and blocking the same shadow cross-section anywhere on it blocks the same amount of irradiance.

**The bug.** The code applied cos(θ_S) twice:

1. **Line 133** of `analysis_Shaded_Irradiance`: `projected_shade_areas` was multiplied by `max(0, sum(-vectors.unit_normal.S .* vectors.unit.p.S, 1))` — i.e., cos(θ_S). The comment said this "gives the shaded areas on the sun," but this projects the cross-section *onto* the surface normal, giving the perpendicular component rather than converting it to actual surface area.

2. **Lines 169/177**: `limb_darkened_projections_SPxS` was formed as `projections_SPxS .* power_Series_Of_Matrix(...)`, meaning cos(θ_S) × LD(cos(θ_S)). The LD factor is correct, but the leading `projections_SPxS` (cos(θ_S)) is a second application of the same factor.

The net effect: the shade irradiance was proportional to cos²(θ_S) × LD(cos(θ_S)) when it should be proportional only to LD(cos(θ_S)). As the observer moves from equator to pole, the shadow shifts ~0.57–0.9 solar radii across the disc. The spurious cos²(θ_S) factor reduced the computed shading by ~31% over this shift, regardless of whether any shade elements actually fell off the disc.

**The fix:**
- **Line 133**: Removed the `.*max(0, sum(-vectors.unit_normal.S .* vectors.unit.p.S, 1))` multiplication entirely. `projected_shade_areas` now correctly represents the shadow cross-section at the sun distance, and nothing more.
- **Line 169**: Changed `limb_darkened_projections_SPxS = projections_SPxS .* power_Series_Of_Matrix(...)` to `limb_darkened_projections_SPxS = power_Series_Of_Matrix(...)`. The plain cos(θ_S) is no longer included; only the limb darkening function LD(cos(θ_S)) is applied.

---

### Bug 2 — Asymmetric limb darkening between the two analysis functions
**Responsible for ~3× underestimate of the absolute shading factor magnitude at all latitudes (no effect on the gradient).**

**Background.** Limb darkening was temporarily commented out in `analysis_Irradiance` during gradient investigation (the TEST comments). It was never commented out in `analysis_Shaded_Irradiance`. The shading factor is computed as:

```
shading_factor = analysis_Shaded_Irradiance / analysis_Irradiance
```

With LD active in the numerator but not the denominator, the denominator (base irradiance) was larger than physically correct. This deflated the shading factor by a roughly constant factor (~3×) at every grid point. Because the deflation was uniform across all latitudes, this bug did not show up in the gradient diagnostic — it was invisible to the tests Morgan ran last night when attempting to isolate the limb darkening effect.

**The fix:** Restored `analysis_Irradiance` to its pre-test state: uncommenting the `limb_darkened_projections_SPxS` calculation and restoring its use in the spectral loop. Both functions now apply limb darkening identically to their respective star-surface calculations.

---

## Test Results Confirming the Fix

The following four variants were computed directly from intermediate values (not through `analysis_General`) to isolate each bug's contribution:

| Variant | Description | Gradient (equator → pole) | SF at equator |
|---------|-------------|--------------------------|---------------|
| A | Current code (both bugs) | **38.40%** | 0.00296 |
| B | Bug 1 fixed only | 7.32% | 0.00300 |
| C | Bug 2 fixed only | 38.40% | 0.00918 |
| D | Both bugs fixed | **7.32%** | **0.00931** |

Key findings:
- Bug 1 alone explains 31 of the 38 percentage points of false gradient.
- Bug 2 has zero effect on the *gradient* but changes the absolute shading factor by ~3×.
- Fixing both (Variant D) gives a 7.32% gradient, consistent with the expected ~7–9% from limb darkening alone.
- The base irradiance (no shade) was confirmed flat to machine precision across the disc partition (Test 01), ruling out any issue in the irradiance-only path.

---

## Secondary Issue (Also Fixed)

In `propagate_Vectors_To_Sphere___M2C___Sr.m`, line 90, the outward unit normal on the sun surface was computed as:
```matlab
vectors.unit_normal.S = t_1 .* vectors.p.H - center_S;
```
The correct expression is:
```matlab
vectors.unit_normal.S = t_1 .* vectors.p.H + vectors.p.s;   % vectors.p.s = vector_P - center_S
```
The intersection point in the common coordinate system is `P + t_1*(H - P)`. The outward normal is `(intersection - center_S) = t_1*(H-P) + (P - center_S) = t_1*vectors.p.H + vectors.p.s`. The previous code omitted the `+ vector_P` offset (absorbed in `vectors.p.s`), producing an angular error of ~0.06–0.14° in the sun surface normal depending on latitude. This was too small to cause the gradient bug, but has been corrected for physical accuracy.
