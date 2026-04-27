# Check Spectrometer Routine — Implementation Specification

> **Purpose of this document**
> This is a complete, software-agnostic specification of the "Check Spectrometer" routine. It is written so that a developer (or an LLM) can re-implement the routine in any environment — Python, MATLAB, LabVIEW, C#, Igor, etc. — without needing to read the original Blick source code.
>
> The routine takes a single measurement of a spectral calibration lamp (typically a mercury vapor lamp), finds the brightest emission line in the spectrum, and fits a modified Gaussian to it. The fit center, width, shape, and quality are reported numerically and visually.

---

## 1. High-Level Overview

### 1.1 What problem does this solve?

Spectrometers drift over time. The pixel-to-wavelength mapping, the optical resolution, and the line-shape function can all change due to thermal expansion, mechanical stress, or detector aging. A spectral lamp emits light at known, narrow, atomic transition wavelengths (e.g. mercury at 253.65 nm, 296.73 nm, 365.02 nm, 404.66 nm, 435.83 nm, 546.07 nm, etc.), so it provides a fixed reference. By repeatedly measuring this lamp and fitting its strongest peak, the user can:

1. Verify the spectrometer is functioning.
2. Track the **center pixel** of a known line over time (wavelength calibration drift).
3. Track the **resolution** (FWHM-like width) over time (focus / slit drift).
4. Track the **line shape exponent** (Gaussian → super-Gaussian → boxy) over time.
5. Catch hardware faults early (saturated pixels, dead pixels, dark-current spikes).

### 1.2 Inputs

| Input | Source | Notes |
|---|---|---|
| Lamp spectrum | Spectrometer hardware | One mean + one std-dev array, one value per detector pixel |
| Spectrometer descriptor | Configuration | Includes total pixel count and number of "blind" (dark) pixels at the end |
| Integration duration | User | Default 3 seconds |
| Saturation target | Constant | 60 % of full scale |
| Half-window size | Constant | `dpix = 40` pixels on each side of the peak |
| Initial resolution guess | Constant | `resol = 5` pixels |
| Fit method ID | Constant | `meth = 11` (modified Gaussian) |

### 1.3 Outputs

| Output | Type | Meaning |
|---|---|---|
| `xcen` | float | Sub-pixel center of the fitted line, in absolute pixel units |
| `resolfit` | float | Fitted width (a Gaussian-σ-like parameter) in pixels |
| `n` (= `a[3]`) | float | Shape exponent of the modified Gaussian (2 = pure Gaussian, >2 = flatter top, <2 = sharper) |
| `rms` | float | Root-mean-square residual of the fit (smaller = better) |
| `err` | int | 0 = converged, 1 = max function evaluations exceeded, 2 = max iterations exceeded |
| Diagnostic plot | image | JPG + native data file, saved to disk |

---

## 2. The Modified Gaussian Model

The fit function used is a **generalized (super-)Gaussian** with a constant background:

```
                                          ⎛  | x - xcen |^n ⎞
    f(x) = B  +  A · exp ⎜ - 0.5 · ⎜ ─────────── ⎟    ⎟
                                          ⎝  |     w     |    ⎠
```

Where:
- `B` — background offset (estimated from the wings of the window).
- `A` — peak amplitude above background.
- `xcen` — sub-pixel line center.
- `w` — width (analogous to σ for `n=2`).
- `n` — shape exponent. `n=2` is a true Gaussian. `n=4`, `n=6`, etc. produce flatter-topped, sharper-edged shapes typical of slit-limited spectrometers.

The fit returns `xcen`, `w`, `n`, plus diagnostic arrays. Typical values for a healthy spectrometer measuring a mercury line:
- `n` between 2 and 6
- `w` between 1 and 8 pixels
- `rms / A` < 1 %

---

## 3. Step-by-Step Procedure

### Step 1 — Adjust Integration Time (auto-exposure)

Before the real measurement, run an auto-exposure loop:

1. Take a short test exposure.
2. Find the maximum signal (excluding any blind/dark pixels).
3. Compute `saturation_fraction = max_signal / detector_full_scale`.
4. Adjust `integration_time_new = integration_time_current × (0.60 / saturation_fraction)`.
5. Repeat until `saturation_fraction` lies in a comfortable band around 60 % (typically 50 – 70 %).

**Important:** during this auto-exposure, do *not* perform dark-current subtraction. The point is to keep the *raw* signal below saturation.

### Step 2 — Acquire the Spectrum

Set:
- Integration time: 3 seconds (overrides Step 1's result for the final acquisition).
- Number of repetitions: 1.

Acquire one measurement. Compute and store, per pixel:
- `signal_mean[i]` — the mean digital count (since `nrep = 1`, this equals the single reading).
- `signal_std[i]` — the standard error / uncertainty estimate.

If any pixel reads at the detector's full-scale value, raise a non-fatal warning ("SATURATION") but continue.

### Step 3 — Dark Correction Using Blind Pixels

Many spectrometers have a small number of physically masked "blind" pixels at the end of the detector array, used to estimate dark current and bias.

The spectrometer descriptor provides two numbers:
- `n_total` — total number of pixels (e.g. 2048).
- `n_blind` — number of blind pixels at the end (e.g. 14). May be zero.

**If `n_blind > 0`:**

```
dark_level = max( signal_mean[ n_total - n_blind  :  n_total ] )
signal_corrected = signal_mean[ 0  :  n_total - n_blind ]  -  dark_level
uncert_corrected = signal_std [ 0  :  n_total - n_blind ]
```

**If `n_blind == 0`:**

```
signal_corrected = signal_mean
uncert_corrected = signal_std
```

From here on, "the spectrum" refers to `signal_corrected` and its length is `N = len(signal_corrected)`.

### Step 4 — Locate the Brightest Peak and Choose a Window

```
indm  = argmax( signal_corrected )            # pixel index of global maximum
dpix  = 40                                    # half-window
ind1  = max( 0, indm - dpix )
ind2  = min( N, indm + dpix + 1 )             # exclusive upper bound
```

Build the windowed arrays:

```
xi    = [ind1, ind1+1, ..., ind2-1]           # absolute pixel indices
yi    = signal_corrected[ind1 : ind2]
uyi   = uncert_corrected[ind1 : ind2]
yimax = signal_corrected[indm]
```

If the user later observes that 40 pixels is too wide (catches a neighboring line) or too narrow (misses the wings), `dpix` is the single tuning knob to adjust.

### Step 5 — Fit the Modified Gaussian

Call a non-linear least-squares fitter with:
- Independent variable: `xi - indm` (centered on the peak so the fit center starts near 0).
- Dependent variable: `yi`.
- Optional weights: `1 / uyi²`.
- Initial guesses: amplitude ≈ `yimax - min(yi)`, center ≈ 0, width ≈ `resol = 5`, exponent ≈ 2, background ≈ `min(yi)`.
- Model: see Section 2.

The fitter must additionally return four classification arrays describing how each window pixel was used:

| Return | Meaning |
|---|---|
| `xxi[0], yyi[0]` | Pixels inside the window but **excluded from the fit** (e.g. the very saturated core, if rejected). May be empty. |
| `xxi[1], yyi[1]` | Pixels classified as **background** (used to estimate `B`). Typically the outermost pixels of the window. |
| `xxi[2], yyi[2]` | Pixels actually **used in the fit** (the body and shoulders of the peak). |
| `xxi[3], yyi[3]` | A **densely sampled evaluation of the fitted curve** for plotting (e.g. 200 points spanning the window). |

The x-coordinates in `xxi[*]` are **relative to `indm`** (so `xxi[3]` runs roughly from `-40` to `+40`). Add `indm` back when plotting against absolute pixel number.

The fitter also returns:
- `err` ∈ {0, 1, 2}
- `a` — full parameter vector. Convention used here: `a[0] = center offset from indm`, `a[3] = shape exponent n`.
- `rms` — RMS residual.
- `resolfit` — fitted width `w`.
- `fitstat`, `fmss` — extra diagnostics, not used downstream.

If your environment has no such fitter, implement Levenberg–Marquardt against the model in Section 2; identify background as the lowest 10 % of the window, and use the rest for fitting.

### Step 6 — Compute Display Quantities

```
xcen   = indm + a[0]                          # absolute sub-pixel center
indbg  = xxi[1] + indm                        # absolute background pixel indices
bg     = mean( signal_corrected[ indbg ] )    # background level for normalization
norm   = signal_corrected - bg
norm   = norm / norm[indm]                    # so the peak == 1.0
pix    = [0, 1, 2, ..., N-1]
```

Build the legend label for the fit:

```
label_fit = "CEN=" + format(xcen, ".2f") +
            ", w=" + format(resolfit, ".2f") +
            ", n=" + format(a[3], ".2f") +
            ", rms=" + format(rms, ".4f")
if err == 1: label_fit += ", >maxfun!"
if err == 2: label_fit += ", >maxiter!"
```

---

## 4. Plot Specification

The plot shows **five overlaid series** on a single Cartesian axis. X-axis is pixel number, Y-axis is normalized signal.

### 4.1 Series Definitions (in draw order)

| # | Name | X data | Y data | Style | Legend text |
|---|---|---|---|---|---|
| 1 | Fitted curve | `xxi[3] + indm` | `yyi[3]` | **Red** solid line, thin (≈0.003 line-weight unit) | "CEN=…, w=…, n=…, rms=…" (see §3.6) |
| 2 | Full normalized spectrum | `pix` (0 … N-1) | `norm` | **Gray** solid line, very thin (≈0.0001) | "SIGNAL" |
| 3 | Data inside window but not fitted | `xxi[0] + indm` | `yyi[0]` | **Light blue** dots, marker size 12, no connecting line | "DATA NOT USED FOR FITTING" |
| 4 | Data actually fitted | `xxi[2] + indm` | `yyi[2]` | **Blue** dots, marker size 1, no connecting line | "FITTING DATA" |
| 5 | Background pixels | `xxi[1] + indm` | `yyi[1]` | **Black** dots, marker size ≈ default ("-1" code), no line | "BACKGROUND DATA" |

> Note on the original code's style codes: the last element of each `lineparsi` is a *style index* into a palette where `0` = red, `-2` = gray, `12` = light blue (large dot), `1` = blue (small dot), `-1` = black (default dot). The first element is the line thickness; `0` means "no connecting line, only markers." When re-implementing, just match the colors and marker/line behavior in the table above.

### 4.2 Y-axis Normalization

All five series share the same Y-axis. Series 1, 3, 4, 5 come straight from the fit/measurement and are **already on a scale where the peak ≈ 1** because the fitter normalizes internally (the modified Gaussian's amplitude `A` is fit against the mean-subtracted, peak-divided window). Series 2 is normalized explicitly in Step 6. The result is a unit-height peak with sensible wings.

### 4.3 Axes

- X-axis label: `"PIXEL"`
- Y-axis label: `"NORMALIZED SIGNAL"`
- X-axis range: see §4.5 (two-pass strategy).
- Y-axis range: see §4.5.
- Linear scale on both axes.

### 4.4 Title

Two lines:

```
{instrument_name} at {location}, routine {routine_name}, {timestamp_human}
Fitting strongest line
```

Example:

```
Pandora123s1 at Innsbruck, routine CheckSpectrometer, 2026-04-27 14:32:11 UTC
Fitting strongest line
```

### 4.5 Two-Pass Axis Limits ("zoom-out reset" + "default narrow view")

The user can press a "reset zoom" key in the GUI; that key returns to the **wide** view. The **default** view shown when the figure first appears is **narrow**. Implement this as two preparation passes:

**Pass 1 — wide limits (saved as the "home" zoom):**
- Compute X limits over **all five series** combined, with a small outward margin (≈1 % each side).
- Compute Y limits over **all five series** combined, with ≈5 % margin.
- Store these as the figure's "initial" / "home" limits.

**Pass 2 — narrow limits (used for display):**
- Compute X limits over series **3, 4, 5 only** (i.e. exclude the full-spectrum gray line and the dense fitted-curve trace; keep only the data points inside the window).
- Compute Y limits over series 3, 4, 5 only.
- Apply these as the visible limits.

The result: by default the user sees a tight zoom on the peak; when they press "reset zoom" they see the entire spectrum with the peak in context.

> Margin formula used by the original `OptimizeAxisLimits`:
> `lo' = lo - margin% × (hi - lo)` and `hi' = hi + margin% × (hi - lo)`.

### 4.6 Legend

- Position: `(x=0.5, y=0.37)` in axes-fraction coordinates — i.e. centered horizontally, in the lower third vertically. (Tweak if the peak in your data crowds this region.)
- Spacing: small inter-row pad (`0.02` in axes-fraction units).
- Frame: hidden (`legend[3] = False`).
- Entries: the five `legtext` strings, in the order above.

### 4.7 Saving

After drawing, save two files to the configured output directory:

```
{instrument_id}{instrument_index}_{location}_{timestamp_compact}_CheckSpectrometer.jpg
{instrument_id}{instrument_index}_{location}_{timestamp_compact}_CheckSpectrometer.ddf
```

- `.jpg` — rasterized figure for human review.
- `.ddf` — native binary/text dump of the underlying X/Y arrays so the plot can be re-rendered later without re-measuring.

`timestamp_compact` is typically `YYYYMMDDThhmmssZ` or similar; pick a format that sorts lexicographically.

---

## 5. Pseudocode (Reference Implementation)

```pseudo
function CheckSpectrometer(spectrometer, location, instrument_name):

    # --- 1. Auto-expose to ~60 % saturation -----------------------
    AdjustIntegrationTime(target = 0.60, dark_correction = false)

    # --- 2. Take the real measurement -----------------------------
    SetSpectrometer(duration_s = 3, n_repetitions = 1)
    raw_mean, raw_std, sat_flag = Measure()
    if sat_flag: Warn("saturation detected")

    # --- 3. Dark correction with blind pixels ---------------------
    n_total, n_blind = spectrometer.pixel_layout
    if n_blind > 0:
        dark = max(raw_mean[n_total - n_blind : n_total])
        signal = raw_mean[0 : n_total - n_blind] - dark
        uncert = raw_std [0 : n_total - n_blind]
    else:
        signal, uncert = raw_mean, raw_std

    # --- 4. Window around the brightest pixel ---------------------
    dpix = 40
    indm = argmax(signal)
    ind1 = max(0, indm - dpix)
    ind2 = min(len(signal), indm + dpix + 1)
    xi   = range(ind1, ind2)
    yi   = signal[ind1 : ind2]
    uyi  = uncert[ind1 : ind2]

    # --- 5. Fit modified Gaussian ---------------------------------
    err, a, rms, resolfit, _, xxi, yyi, _ = FitPeak(
        xi, yi, uyi,
        center_pixel = indm,
        initial_width = 5,
        method = "modified_gaussian"
    )

    # --- 6. Normalize full spectrum for display -------------------
    indbg = xxi[1] + indm
    bg    = mean(signal[indbg])
    norm  = (signal - bg) / (signal[indm] - bg)
    pix   = range(0, len(signal))

    # --- 7. Build five plot series --------------------------------
    xcen = indm + a[0]
    label_fit = format("CEN=%.2f, w=%.2f, n=%.2f, rms=%.4f",
                       xcen, resolfit, a[3], rms)
    if err == 1: label_fit += ", >maxfun!"
    if err == 2: label_fit += ", >maxiter!"

    series = [
      Series(x = xxi[3] + indm, y = yyi[3], style = "red line thin",
             label = label_fit),
      Series(x = pix,           y = norm,   style = "gray line very thin",
             label = "SIGNAL"),
      Series(x = xxi[0] + indm, y = yyi[0], style = "lightblue dots big",
             label = "DATA NOT USED FOR FITTING"),
      Series(x = xxi[2] + indm, y = yyi[2], style = "blue dots small",
             label = "FITTING DATA"),
      Series(x = xxi[1] + indm, y = yyi[1], style = "black dots default",
             label = "BACKGROUND DATA"),
    ]

    # --- 8. Two-pass axes -----------------------------------------
    title = instrument_name + " at " + location + ", routine "
          + "CheckSpectrometer, " + Now() + "\nFitting strongest line"

    wide_x = AxisLimits(all_series,  margin = 1 %)
    wide_y = AxisLimits(all_series,  margin = 5 %)
    DrawFigure(series, xlim = wide_x, ylim = wide_y,
               xlabel = "PIXEL", ylabel = "NORMALIZED SIGNAL",
               title = title, save_as_home_zoom = true)

    narrow_x = AxisLimits(series[2:], margin = default)
    narrow_y = AxisLimits(series[2:], margin = default)
    DrawFigure(series, xlim = narrow_x, ylim = narrow_y,
               xlabel = "PIXEL", ylabel = "NORMALIZED SIGNAL",
               title = title, save_as_home_zoom = false)

    # --- 9. Persist ----------------------------------------------
    base = OutputDir() + InstrumentTag() + "_" + location + "_"
         + Timestamp("compact") + "_CheckSpectrometer"
    Save(base + ".jpg")
    Save(base + ".ddf")

    return { xcen, resolfit, n: a[3], rms, err }
```

---

## 6. Designing It as a Service

If you wrap this routine as a network/REST service, here is a clean interface:

### 6.1 Endpoint

```
POST /spectrometer/check
```

### 6.2 Request body

```json
{
  "duration_s": 3,
  "saturation_target": 0.60,
  "half_window_pixels": 40,
  "initial_width_pixels": 5,
  "fit_method": "modified_gaussian",
  "instrument_id": "Pandora123s1",
  "location": "Innsbruck"
}
```

### 6.3 Response body

```json
{
  "status": "ok",
  "fit": {
    "center_pixel":   1024.37,
    "width_pixels":   3.12,
    "shape_exponent": 2.84,
    "rms_residual":   0.0017,
    "converged":      true,
    "warnings":       []
  },
  "saturation_detected": false,
  "dark_level_counts":   123.4,
  "plot": {
    "jpg_path": "/data/2026/04/27/Pandora123s1_Innsbruck_20260427T143211Z_CheckSpectrometer.jpg",
    "ddf_path": "/data/2026/04/27/Pandora123s1_Innsbruck_20260427T143211Z_CheckSpectrometer.ddf"
  },
  "data": {
    "pixel":           [...],
    "signal_norm":     [...],
    "fit_curve_pixel": [...],
    "fit_curve_value": [...]
  }
}
```

### 6.4 Error semantics

| Condition | HTTP | `status` | Notes |
|---|---|---|---|
| Spectrometer unreachable | 503 | `"hardware_error"` | retry once, then fail |
| Saturation persists after auto-expose | 200 | `"ok"` | with warning `"saturation"` |
| Fit did not converge | 200 | `"ok"` | `converged: false`, still return numbers |
| Peak amplitude below noise floor | 422 | `"no_peak"` | refuse to fit |
| Invalid request | 400 | `"bad_request"` | |

### 6.5 Concurrency

The spectrometer is a single-user resource — serialize calls behind a mutex. Each call typically takes `auto_expose_time + duration_s + fit_time ≈ 5 – 15 s`.

---

## 7. Validation Tests

A re-implementation can be considered correct if, given the same input arrays, it produces:

1. Same `indm` (exact integer match).
2. Same dark-corrected spectrum (exact float match modulo rounding).
3. `xcen`, `resolfit`, `n` within 0.1 % of the reference implementation.
4. `rms` within 1 % of the reference implementation.
5. Plot containing all five series with the colors, markers, and labels of §4.1.
6. Default view zoomed to series 3 – 5; "home" view spanning all series.
7. JPG + DDF files saved with the correct naming convention.

A canned mercury-lamp dataset and the reference outputs should be checked into the test suite.

---

## 8. Glossary

| Term | Meaning |
|---|---|
| **Blind pixels** | Detector pixels physically masked from light, used to estimate dark current. |
| **Dark correction** | Subtracting the dark/bias level from each pixel. |
| **Modified Gaussian** | Generalized Gaussian `exp(-0.5·|Δx/w|^n)` with variable exponent `n`. |
| **Resolution (`w`)** | Width parameter of the fitted line, in pixels. Convertible to nm via the wavelength-per-pixel dispersion. |
| **DDF** | Native data dump file format used by the original software; treat as "structured array dump." |
| **`indm`** | Integer pixel index of the maximum signal. |
| **`xcen`** | Sub-pixel center of the fitted line, in absolute pixel units. |
