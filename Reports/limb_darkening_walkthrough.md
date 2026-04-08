# Limb Darkening in the PSF Sunshade Model: Step-by-Step Walkthrough

**Planetary Sunshade Foundation — Analysis Document**
**March 20, 2026**

---

## 1. What limb darkening IS (physical intuition)

The Sun is not uniformly bright. When you look at the center of the solar disk, your line of sight penetrates deep into the hot, dense photosphere. When you look near the edge (the "limb"), your line of sight skims through the upper, cooler layers. Since radiation intensity scales with T⁴ (Stefan-Boltzmann), the cooler gas at the limb emits substantially less light. The net effect: the limb appears dimmer than the center — about 61% as bright for total (bolometric) intensity, and as low as 30% at visible wavelengths like 550 nm.

For a sunshade constellation, this matters because a shade element blocking part of the solar disk near the center blocks MORE flux per unit area than the same shade element blocking near the limb.

---

## 2. The mathematical model (PDF Section 1.4, Equation 1.4.1)

The PSF model expresses the observed intensity at angle θ from the surface normal as a power series in cos θ:

```
I(θ) = I(0) × Σ aᵢ cosⁱ(θ)       [Eq. 1.4.1]
```

where:
- θ is the angle of incidence at the star's surface (0° at disk center, 90° at the limb)
- I(0) is the central intensity
- The aᵢ coefficients must sum to 1 (so that I(0)/I(0) = 1 at center)
- The sum starts from i = 0 (a constant term is permitted)

This is more general than the common "linear limb darkening law" I(μ) = I(0)[1 − u(1 − μ)] found on Wikipedia. The two are related: the Wikipedia form with u = 0.39 is equivalent to a₀ = 1 − u = 0.61, a₁ = u = 0.39, which is exactly what the PSF model uses.

---

## 3. The actual coefficients from the Excel file

From `star spectral data.xlsx`, the default (total intensity) band uses four coefficient columns:

| Column | Coefficient | Value | Maps to Eq. 1.4.1 as | Role in I(θ)/I(0) |
|--------|-------------|-------|-----------------------|---------------------|
| C1     | a₀          | 0.61  | constant term         | 0.61 × 1            |
| C2     | a₁          | 0.39  | cos θ term            | 0.39 × cos θ        |
| C3     | a₂          | 0.00  | cos²θ term            | 0 × cos²θ           |
| C4     | a₃          | 0.00  | cos³θ term            | 0 × cos³θ           |

These are cited from Tripathi et al. (2020) for total solar intensity.

So the limb darkening profile is:

```
I(θ) / I(0) = 0.61 + 0.39 · cos(θ)
```

**Key values:**
- At disk center (θ = 0°): I/I(0) = 0.61 + 0.39 = **1.00** ✓
- At the limb (θ = 90°): I/I(0) = 0.61 + 0 = **0.61** (limb is 61% as bright as center)

This is physically reasonable for total (all-wavelength) solar intensity. At individual visible wavelengths, the darkening is stronger — e.g. at 550 nm the limb drops to about 30% of center.

---

## 4. How the code maps these coefficients through `power_Series_Of_Matrix`

The function `power_Series_Of_Matrix` (in `operations/`) computes a polynomial from an input matrix and a vector of coefficients. Tracing through the loop:

```
i=1:  matrix_power = ones(...)     → cos⁰θ = 1     → adds C1 × 1
i=2:  matrix_power *= matrix       → cos¹θ = cosθ   → adds C2 × cosθ
i=3:  matrix_power *= matrix       → cos²θ           → adds C3 × cos²θ
i=4:  matrix_power *= matrix       → cos³θ           → adds C4 × cos³θ
```

So `power_series(cos θ, [C1, C2, C3, C4])` computes:

```
C1 + C2·cosθ + C3·cos²θ + C4·cos³θ
```

**This is exactly Equation 1.4.1's** Σ aᵢ cosⁱ(θ), with C1 = a₀, C2 = a₁, C3 = a₂, C4 = a₃.

For the default coefficients: `0.61 + 0.39·cosθ`

---

## 5. How the code normalizes these coefficients

### 5a. Why normalization is needed

The raw coefficients describe the *shape* of the angular intensity profile. But the code needs to use them with a known total luminosity (L☉ = 3.828 × 10²⁶ W). The normalization ensures that integrating the intensity over all emission directions from a surface patch recovers the correct total power.

### 5b. The normalization integral (`import_Spectral_Data`, lines 68-72)

The code computes:

```matlab
spectral_data.normalizing_coefficients = 8*C4 + 3*pi*C3 + 12*C2 + 6*pi*C1;
spectral_data.normalizing_coefficients = spectral_data.normalizing_coefficients / 6;
```

Each multiplier comes from the **Wallis integral** — the analytic result of integrating cosⁿθ over [−π/2, π/2]:

| Power series term | Integral ∫₋π/₂^{π/2} cosⁿθ dθ | × 6 (code form) |
|-------------------|-------------------------------|------------------|
| C1 × cos⁰θ = C1 × 1   | C1 × **π**     | **6π** · C1  |
| C2 × cos¹θ             | C2 × **2**     | **12** · C2  |
| C3 × cos²θ             | C3 × **π/2**   | **3π** · C3  |
| C4 × cos³θ             | C4 × **4/3**   | **8**  · C4  |

### 5c. How the Wallis integrals arise

Each integral is derived from:

```
∫₋π/₂^{π/2} cosⁿθ dθ
```

- **n = 0:** integrating 1 over [−π/2, π/2] gives **π**
- **n = 1:** ∫cos θ dθ = sin θ |₋π/₂^{π/2} = 1 − (−1) = **2**
- **n = 2:** using cos²θ = (1 + cos 2θ)/2 and integrating gives **π/2**
- **n = 3:** using cos³θ = cos θ(1 − sin²θ) and substitution gives **4/3**

The factor of 6 in the denominator is a common scaling factor so the code can work with integers and simple multiples of π.

### 5d. The result for the default coefficients

```
N = (6π × 0.61 + 12 × 0.39 + 0 + 0) / 6
  = (11.498 + 4.68) / 6
  = 2.6964
```

### 5e. The normalized coefficients

```
C1_norm = 0.61 / 2.6964 = 0.22623
C2_norm = 0.39 / 2.6964 = 0.14464
```

After normalization, the coefficients no longer sum to 1 — they sum to 0.37087. This is correct: the normalization absorbs the relationship between "intensity profile shape" and "total luminosity" so that the downstream computation produces correct W/m² values. The limb/center ratio is preserved: 0.22623 / 0.37087 = 0.61, as expected.

---

## 6. How the code applies limb darkening to irradiance

### 6a. The irradiance calculation (`analysis_Irradiance`, lines 60 and 101)

For each star surface patch (Sα) and planet surface patch (Sβ), the computation proceeds in two stages:

**Line 60 — everything except limb darkening and star-surface projection:**

```
intermediate = L☉ × area_ratio × (1/r²) × cos(planet angle)
```

**Line 101 — limb darkening and star-surface geometric projection:**

```matlab
limb_darkened = projections_SPxS .* power_Series_Of_Matrix(projections_SPxS, coefficients);
```

Here `projections_SPxS` = cos θ (the dot product of the star surface normal with the star-to-planet direction). This line multiplies two things together:

1. **cos θ** — the geometric projection of the star surface element (Lambert's cosine law for foreshortening)
2. **power_series(cos θ)** = C1_n + C2_n · cos θ — the limb darkening profile

The full angular factor becomes:

```
cos θ × [C1_n + C2_n · cos θ] = C1_n · cos θ + C2_n · cos²θ
```

**The first cos θ is NOT part of limb darkening** — it is the geometric foreshortening of the star surface patch when viewed off-axis. Even a uniformly bright star has this factor. The power series alone is the limb darkening.

### 6b. Line 109 — assembling the result

```matlab
spectral_irradiance = Iν × limb_darkened × intermediate;
```

This multiplies the spectral band fraction (Iν), the angular factors (limb-darkened projection), and the intermediate result (luminosity, area ratio, 1/r², planet projection).

---

## 7. How limb darkening applies to shaded irradiance

### 7a. The shaded irradiance calculation (`analysis_Shaded_Irradiance`, lines 152-165)

For each planet-shade pair, the code:

1. Projects the shade onto the star surface (finding where it falls, at some angle θ)
2. Computes `projections_SPxS` — the cosine of the incidence angle at that star-surface point
3. Applies the identical limb darkening formula:

```matlab
limb_darkened_projections_SPxS = projections_SPxS .* power_Series_Of_Matrix(projections_SPxS, coefficients);
```

This means the flux blocked by a shade element depends on WHERE on the solar disk it falls. The limb darkening profile I(θ)/I(0) = 0.61 + 0.39·cos θ means:

- Shade at disk center (θ = 0°): blocks flux proportional to intensity **1.00**
- Shade at r/R = 0.5 (θ = 30°): blocks flux proportional to intensity **0.95**
- Shade at r/R = 0.87 (θ = 60°): blocks flux proportional to intensity **0.81**
- Shade at the limb (θ = 90°): blocks flux proportional to intensity **0.61**

### 7b. Quantitative comparison: shade effectiveness with vs. without limb darkening

The table below shows the relative flux blocked by a shade at each position, comparing the limb-darkened model to a uniform (Lambertian) baseline. The "Blocked flux" columns show the flux blocked as a percentage of the maximum (shade at disk center), incorporating both the intensity profile and the geometric cos θ projection.

| θ (deg) | r/R  | cos θ | I(θ) with LD | Blocked flux (LD) | Blocked flux (no LD) | LD / no-LD |
|---------|------|-------|-------------|-------------------|---------------------|------------|
| 0°      | 0.00 | 1.000 | 1.000       | 100.0%            | 100.0%              | 1.165      |
| 10°     | 0.17 | 0.985 | 0.994       | 97.9%             | 97.0%               | 1.158      |
| 20°     | 0.34 | 0.940 | 0.977       | 91.8%             | 88.3%               | 1.138      |
| 30°     | 0.50 | 0.866 | 0.948       | 82.1%             | 75.0%               | 1.104      |
| 45°     | 0.71 | 0.707 | 0.886       | 62.6%             | 50.0%               | 1.032      |
| 60°     | 0.87 | 0.500 | 0.805       | 40.3%             | 25.0%               | 0.938      |
| 70°     | 0.94 | 0.342 | 0.743       | 25.4%             | 11.7%               | 0.866      |
| 80°     | 0.98 | 0.174 | 0.678       | 11.8%             | 3.0%                | 0.790      |
| 85°     | 1.00 | 0.087 | 0.644       | 5.6%              | 0.8%                | 0.750      |

**The crossover point is near θ ≈ 56° (r/R ≈ 0.83).** Inside this radius, limb darkening makes the shade MORE effective than uniform-disk math predicts. Outside, it makes the shade LESS effective.

---

## 8. What "no limb darkening" means in this code

To run without limb darkening, set in the Excel file:

```
C1 = 1, C2 = 0, C3 = 0, C4 = 0
```

This gives I(θ)/I(0) = 1 (uniform brightness in all directions). The power series returns a constant 1, and the geometric cos θ projection on line 101 provides the standard Lambertian emission pattern. The normalization would be N = π.

---

## 9. Summary: where counterintuitive results may arise

### 9a. The LD / no-LD crossover

With limb darkening, shading near the center (r/R < 0.83) blocks MORE flux than uniform-disk math predicts. Shading near the edge (r/R > 0.83) blocks LESS. If your constellation is concentrated near the optical axis (as most practical designs are), limb darkening will make the shade slightly more effective overall.

### 9b. The constant term dominates

Because a₀ = 0.61 is the dominant coefficient, the intensity profile is relatively flat — varying only from 1.00 at center to 0.61 at the limb. The geometric cos θ projection factor actually has a stronger effect on shade efficiency than limb darkening does.

### 9c. Total irradiance is slightly higher with LD

Integrating the full expression cos θ × (0.61 + 0.39·cos θ) over the hemisphere gives a slightly higher value (~5%) than the no-LD case cos θ × (1/π), when both are normalized to the same total luminosity. This is a real physical effect: limb darkening concentrates emission toward the normal direction, which is more aligned with the planet.

### 9d. Normalization is critical

If the coefficients were used without the normalization step (dividing by N = 2.6964), all irradiance values would be scaled incorrectly.
