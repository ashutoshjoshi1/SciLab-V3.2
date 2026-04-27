"""Check Spectrometer service.

Shines a spectral lamp into the spectrometer, records its output, finds the
brightest peak, fits a modified Gaussian (generalized bell-curve), and
produces diagnostic plots + CSV results.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

LOGGER = logging.getLogger(__name__)

# Half-window around the peak in pixels
_HALF_WIN = 40
# Auto-IT target fraction of detector max (60%)
_AUTO_IT_TARGET_FRAC = 0.60
# Maximum auto-IT iterations
_AUTO_IT_MAX_ITERS = 50


# ---------------------------------------------------------------------------
# Modified Gaussian (generalized bell-curve)
#   f(x) = exp( -0.5 * |x / a2|^(2*a3) )
# Parameters: [center_offset, amplitude, a2 (width), a3 (shape exponent)]
# ---------------------------------------------------------------------------

def _mgauss(x: np.ndarray, c: float, amp: float, a2: float, a3: float) -> np.ndarray:
    """Modified Gaussian with variable shape exponent."""
    a2 = max(abs(a2), 1e-6)
    a3 = max(abs(a3), 0.1)
    return amp * np.exp(-0.5 * np.abs((x - c) / a2) ** (2.0 * a3))


def _fit_mgauss(
    x: np.ndarray, y: np.ndarray, sigma: Optional[np.ndarray] = None
) -> Tuple[int, np.ndarray, float, float, np.ndarray, np.ndarray]:
    """
    Fit a modified Gaussian to (x, y).

    Returns
    -------
    err       : 0=success, 1=max_fev, 2=max_iter
    a         : [center, amplitude, a2, a3]
    rms       : root-mean-square residual
    resolfit  : fitted FWHM in pixels (derived from a2, a3)
    """
    try:
        from scipy.optimize import curve_fit  # type: ignore
    except ImportError as exc:
        raise RuntimeError("scipy is required for Check Spectrometer fit") from exc

    if len(x) < 4:
        return 2, np.array([0.0, 1.0, 5.0, 1.0]), np.inf, 10.0

    peak_idx = int(np.argmax(y))
    amp0 = float(y[peak_idx])
    c0 = float(x[peak_idx])
    a2_0 = max(float(np.std(x)) * 0.5, 3.0)

    p0 = [c0, amp0, a2_0, 1.0]
    bounds_lo = [x[0] - 1, 0, 0.5, 0.1]
    bounds_hi = [x[-1] + 1, amp0 * 2, (x[-1] - x[0]) * 0.6, 10.0]

    err = 0
    try:
        popt, _ = curve_fit(
            _mgauss,
            x,
            y,
            p0=p0,
            sigma=sigma,
            bounds=(bounds_lo, bounds_hi),
            maxfev=5000,
        )
        a = popt
    except RuntimeError:
        err = 1
        a = np.array(p0, dtype=float)
    except Exception:
        err = 2
        a = np.array(p0, dtype=float)

    y_fit = _mgauss(x, *a)
    rms = float(np.sqrt(np.mean((y - y_fit) ** 2)))

    # FWHM: solve exp(-0.5 * |dx/a2|^(2*a3)) = 0.5
    # |dx/a2|^(2*a3) = 2*ln(2)  =>  dx = a2 * (2*ln(2))^(1/(2*a3))
    a2, a3 = abs(a[2]), max(abs(a[3]), 0.1)
    resolfit = 2.0 * a2 * (2.0 * np.log(2.0)) ** (1.0 / (2.0 * a3))

    return err, a, rms, resolfit


# ---------------------------------------------------------------------------
# Public result dataclass
# ---------------------------------------------------------------------------

@dataclass
class CheckSpectrometerResult:
    xcen: float
    resolfit: float
    shape_exponent: float
    rms: float
    fit_err: int
    auto_it_ms: float
    plot_path: str
    csv_path: str
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Service
# ---------------------------------------------------------------------------

class CheckSpectrometerService:
    """Run the Check Spectrometer routine on a connected spectrometer."""

    def __init__(self, output_dir: Path, instrument_name: str = "", location: str = ""):
        self.output_dir = Path(output_dir)
        self.instrument_name = instrument_name
        self.location = location

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self, spec) -> CheckSpectrometerResult:
        """
        Execute the full check routine.

        Parameters
        ----------
        spec : SpectrometerBackend
            A connected spectrometer instance with attributes:
              - eff_saturation_limit (float)
              - npix_active (int)
              - rcm (np.ndarray)
              - rcs (np.ndarray)
              - npix_blind_left (int, optional)
              - npix_blind_right (int, optional)
            and methods: set_it(ms), measure(ncy).
        """
        warnings: List[str] = []
        sat_limit = float(getattr(spec, "eff_saturation_limit", 65535))

        # 1. Auto-adjust integration time to ~60% of saturation
        auto_it_ms = self._auto_it(spec, sat_limit, warnings)

        # 2. Measure spectrum (3-second integration, single run)
        spec.set_it(3000.0)
        result = spec.measure(ncy=1)
        if result != "OK":
            LOGGER.warning("Check spectrometer measure returned: %s", result)

        signal_mean: np.ndarray = np.array(spec.rcm, dtype=float)
        signal_std: np.ndarray = np.array(spec.rcs, dtype=float)

        if np.any(signal_mean >= sat_limit):
            warnings.append("WARNING: one or more pixels are saturated.")
            LOGGER.warning("Saturated pixels detected during check spectrometer")

        # 3. Raw data + blind-pixel dark correction
        signal_mean, signal_std = self._apply_blind_correction(
            spec, signal_mean, signal_std
        )

        # 4. Window around brightest peak
        win_data = self._extract_peak_window(signal_mean, signal_std)

        # 5. Fit modified Gaussian
        fit_err, a, rms, resolfit = _fit_mgauss(
            win_data["x_win"].astype(float),
            win_data["y_win"],
            sigma=np.where(win_data["s_win"] > 0, win_data["s_win"], None),
        )

        xcen = win_data["indm"] + a[0]

        if fit_err == 1:
            warnings.append("NOTE: fit hit max function evaluations.")
        if fit_err == 2:
            warnings.append("NOTE: fit hit max iterations.")

        # 6-8. Build plot
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sn = str(getattr(spec, "sn", "Unknown")).strip()
        plot_path = self._build_plot(
            signal_mean,
            win_data,
            a,
            rms,
            resolfit,
            xcen,
            fit_err,
            sn,
            timestamp,
            warnings,
        )

        # 9. Save CSV
        csv_path = self._save_csv(
            signal_mean,
            win_data,
            a,
            rms,
            resolfit,
            xcen,
            auto_it_ms,
            sn,
            timestamp,
        )

        return CheckSpectrometerResult(
            xcen=xcen,
            resolfit=resolfit,
            shape_exponent=float(a[3]),
            rms=rms,
            fit_err=fit_err,
            auto_it_ms=auto_it_ms,
            plot_path=str(plot_path),
            csv_path=str(csv_path),
            warnings=warnings,
        )

    # ------------------------------------------------------------------
    # Step 1: auto-IT
    # ------------------------------------------------------------------

    def _auto_it(self, spec, sat_limit: float, warnings: List[str]) -> float:
        target = _AUTO_IT_TARGET_FRAC * sat_limit
        current_it = float(getattr(spec, "it_ms", 10.0))
        if current_it <= 0:
            current_it = 10.0

        for _ in range(_AUTO_IT_MAX_ITERS):
            spec.set_it(current_it)
            res = spec.measure(ncy=1)
            if res != "OK":
                break

            rcm = np.asarray(spec.rcm, dtype=float)
            if rcm.size == 0:
                break
            peak = float(np.max(rcm))
            if peak <= 0:
                current_it = min(current_it * 2.0, 3000.0)
                continue

            ratio = target / peak
            new_it = current_it * ratio

            # Clamp
            new_it = max(0.2, min(new_it, 3000.0))

            if abs(new_it - current_it) / max(current_it, 1e-6) < 0.05:
                current_it = new_it
                break
            current_it = new_it

        spec.set_it(current_it)
        LOGGER.info("Check spectrometer auto-IT settled at %.2f ms", current_it)
        return current_it

    # ------------------------------------------------------------------
    # Step 3: blind pixel correction
    # ------------------------------------------------------------------

    def _apply_blind_correction(
        self,
        spec,
        signal_mean: np.ndarray,
        signal_std: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        npix_blind_right = int(getattr(spec, "npix_blind_right", 0))
        npix_blind_left = int(getattr(spec, "npix_blind_left", 0))

        if npix_blind_right > 0:
            blind_vals = signal_mean[-npix_blind_right:]
            dark_level = float(np.max(blind_vals))
            LOGGER.info(
                "Dark correction from %d right blind pixels: dark=%.1f",
                npix_blind_right,
                dark_level,
            )
            signal_mean = signal_mean[:-npix_blind_right] - dark_level
            signal_std = signal_std[:-npix_blind_right]
        elif npix_blind_left > 0:
            blind_vals = signal_mean[:npix_blind_left]
            dark_level = float(np.max(blind_vals))
            LOGGER.info(
                "Dark correction from %d left blind pixels: dark=%.1f",
                npix_blind_left,
                dark_level,
            )
            signal_mean = signal_mean[npix_blind_left:] - dark_level
            signal_std = signal_std[npix_blind_left:]

        return signal_mean, signal_std

    # ------------------------------------------------------------------
    # Step 4: peak window extraction
    # ------------------------------------------------------------------

    def _extract_peak_window(
        self, signal_mean: np.ndarray, signal_std: np.ndarray
    ) -> dict:
        npix = len(signal_mean)
        indm = int(np.argmax(signal_mean))
        peak_signal = float(signal_mean[indm])

        lo = max(0, indm - _HALF_WIN)
        hi = min(npix - 1, indm + _HALF_WIN)

        x_win = np.arange(lo, hi + 1)
        y_win = signal_mean[lo : hi + 1]
        s_win = signal_std[lo : hi + 1]

        return {
            "indm": indm,
            "peak_signal": peak_signal,
            "lo": lo,
            "hi": hi,
            "x_win": x_win,
            "y_win": y_win,
            "s_win": s_win,
        }

    # ------------------------------------------------------------------
    # Steps 6-8: build the 5-series diagnostic plot
    # ------------------------------------------------------------------

    def _build_plot(
        self,
        signal_mean: np.ndarray,
        win_data: dict,
        a: np.ndarray,
        rms: float,
        resolfit: float,
        xcen: float,
        fit_err: int,
        sn: str,
        timestamp: str,
        warnings: List[str],
    ) -> Path:
        from matplotlib.figure import Figure

        npix = len(signal_mean)
        indm = win_data["indm"]
        lo = win_data["lo"]
        hi = win_data["hi"]
        x_win = win_data["x_win"]
        y_win = win_data["y_win"]
        peak_signal = win_data["peak_signal"]

        # ----------------------------------------------------------
        # Background estimate: average of outermost 25% on each side
        # ----------------------------------------------------------
        bg_outer = max(1, len(x_win) // 4)
        bg_pixels = np.concatenate([x_win[:bg_outer], x_win[-bg_outer:]])
        bg_vals = np.array([signal_mean[int(p)] for p in bg_pixels if 0 <= int(p) < npix])
        bg_level = float(np.mean(bg_vals)) if len(bg_vals) > 0 else 0.0
        peak_above_bg = peak_signal - bg_level
        norm_denom = peak_above_bg if peak_above_bg > 0 else max(float(np.max(np.abs(y_win))), 1.0)

        # ----------------------------------------------------------
        # Threshold-based split of window into 3 groups
        #   background   : normalised signal < 0.05
        #   not-used     : normalised signal > 0.95  (very core)
        #   fitting      : everything in between
        # ----------------------------------------------------------
        y_win_norm = (y_win - bg_level) / norm_denom
        mask_bg       = y_win_norm < 0.05
        mask_not_used = y_win_norm > 0.95
        mask_fitting  = ~mask_bg & ~mask_not_used

        # ----------------------------------------------------------
        # Full spectrum normalised (Series 2)
        # ----------------------------------------------------------
        signal_sub  = signal_mean - bg_level
        signal_norm = signal_sub / norm_denom
        x_all = np.arange(npix)

        # ----------------------------------------------------------
        # Fitted curve (Series 1) – evaluated over the window
        # ----------------------------------------------------------
        x_fit_local = np.linspace(float(x_win[0]) - indm, float(x_win[-1]) - indm, 500)
        y_fit_raw   = _mgauss(x_fit_local, a[0], a[1], a[2], a[3])
        amp = a[1] if a[1] > 0 else 1.0
        y_fit_norm  = y_fit_raw / amp
        x_fit_abs   = x_fit_local + indm

        # ----------------------------------------------------------
        # Series labels
        # ----------------------------------------------------------
        lbl_warn = "  [iter limit]" if fit_err in (1, 2) else ""
        label1 = f"CEN={xcen:.2f}, w={resolfit:.2f}, n={a[3]:.2f}, rms={rms:.4f}{lbl_warn}"
        label2 = "SIGNAL"
        label3 = "DATA NOT USED FOR FITTING"
        label4 = "FITTING DATA"
        label5 = "BACKGROUND DATA"

        COLOR1 = "red"
        COLOR2 = "#555555"   # dark gray for signal line
        COLOR3 = "lightblue"
        COLOR4 = "blue"
        COLOR5 = "black"

        # ----------------------------------------------------------
        # Build figure – reference style
        # ----------------------------------------------------------
        fig = Figure(figsize=(14, 7))
        ax  = fig.add_subplot(111)

        # Gray background, dashed white grid (matches reference)
        ax.set_facecolor("#c8c8c8")
        fig.patch.set_facecolor("#c8c8c8")
        ax.grid(True, linestyle="--", linewidth=0.7, color="white", alpha=0.9, zorder=0)

        # Series 2 — full normalised signal (thin dark-gray line)
        ax.plot(x_all, signal_norm, color=COLOR2, lw=0.9, label=label2, zorder=1)

        # Series 5 — background dots (black)
        if np.any(mask_bg):
            ax.plot(x_win[mask_bg], y_win_norm[mask_bg],
                    ".", color=COLOR5, ms=6, label=label5, zorder=2)

        # Series 3 — not-used dots (light blue)
        if np.any(mask_not_used):
            ax.plot(x_win[mask_not_used], y_win_norm[mask_not_used],
                    ".", color=COLOR3, ms=7, label=label3, zorder=3)

        # Series 4 — fitting dots (blue)
        if np.any(mask_fitting):
            ax.plot(x_win[mask_fitting], y_win_norm[mask_fitting],
                    ".", color=COLOR4, ms=7, label=label4, zorder=4)

        # Series 1 — fitted curve (red)
        ax.plot(x_fit_abs, y_fit_norm, "-", color=COLOR1, lw=2.2, label=label1, zorder=5)

        # ----------------------------------------------------------
        # Axes labels
        # ----------------------------------------------------------
        ax.set_xlabel("PIXEL", fontsize=13, fontweight="bold")
        ax.set_ylabel("NORMALIZED SIGNAL", fontsize=13, fontweight="bold")

        # ----------------------------------------------------------
        # Two-line title: "instrument at location, routine CS, datetime"
        #                 "Fitting strongest line"
        # ----------------------------------------------------------
        now_str  = datetime.now().strftime("%a %d %b %Y, %H:%M:%S")
        location = getattr(self, "location", "") or ""
        if location:
            line1 = f"{self.instrument_name} at {location}, routine CS, {now_str}"
        else:
            line1 = f"{self.instrument_name}, routine CS, {now_str}"
        ax.set_title(f"{line1}\nFitting strongest line", fontsize=11)

        # ----------------------------------------------------------
        # Legend – upper left, no frame box, text coloured by series
        # ----------------------------------------------------------
        legend_order   = [label1, label2, label3, label4, label5]
        legend_colors  = [COLOR1, COLOR2, COLOR3, COLOR4, COLOR5]
        handles_dict   = {h.get_label(): h for h in ax.get_lines()}

        ordered_handles = []
        ordered_labels  = []
        for lbl in legend_order:
            h = handles_dict.get(lbl)
            if h is not None:
                ordered_handles.append(h)
                ordered_labels.append(lbl)

        leg = ax.legend(
            ordered_handles, ordered_labels,
            loc="upper left",
            fontsize=9,
            frameon=False,
        )
        for text, color in zip(leg.get_texts(), legend_colors[: len(ordered_labels)]):
            text.set_color(color)
            text.set_fontweight("bold")

        # ----------------------------------------------------------
        # Pass 1 – wide x-limits (full detector view)
        # ----------------------------------------------------------
        ax.set_xlim(-5, npix + 5)
        all_y = np.concatenate([signal_norm, y_win_norm, y_fit_norm])
        y_lo  = float(np.nanmin(all_y))
        y_hi  = float(np.nanmax(all_y))
        ax.set_ylim(min(y_lo - 0.05, -0.05), max(y_hi + 0.05, 1.1))

        # ----------------------------------------------------------
        # Pass 2 – narrow default view (zoom around peak, Series 3-5)
        # ----------------------------------------------------------
        pad_x = max(5, int((hi - lo) * 0.20))
        ax.set_xlim(float(x_win[0]) - pad_x, float(x_win[-1]) + pad_x)

        narrow_y = np.concatenate([y_win_norm, y_fit_norm])
        y_lo_n   = float(np.nanmin(narrow_y))
        y_hi_n   = float(np.nanmax(narrow_y))
        ax.set_ylim(min(y_lo_n - 0.05, -0.05), max(y_hi_n + 0.07, 1.1))

        fig.tight_layout()

        # ----------------------------------------------------------
        # Save as JPEG
        # ----------------------------------------------------------
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plot_path = self.output_dir / f"CheckSpectrometer_{sn}_{timestamp}.jpg"
        fig.savefig(plot_path, dpi=150, bbox_inches="tight",
                    format="jpeg", pil_kwargs={"quality": 92})
        LOGGER.info("Check spectrometer plot saved: %s", plot_path)
        return plot_path

    # ------------------------------------------------------------------
    # Step 9: save CSV
    # ------------------------------------------------------------------

    def _save_csv(
        self,
        signal_mean: np.ndarray,
        win_data: dict,
        a: np.ndarray,
        rms: float,
        resolfit: float,
        xcen: float,
        auto_it_ms: float,
        sn: str,
        timestamp: str,
    ) -> Path:
        try:
            import pandas as pd  # type: ignore
        except ImportError as exc:
            raise RuntimeError("pandas is required to save Check Spectrometer CSV") from exc

        npix = len(signal_mean)
        indm = win_data["indm"]
        lo = win_data["lo"]
        hi = win_data["hi"]

        # Summary row
        summary = {
            "SN": sn,
            "Timestamp": timestamp,
            "AutoIT_ms": round(auto_it_ms, 4),
            "MeasureIT_ms": 3000.0,
            "xcen_pixel": round(xcen, 3),
            "resolfit_pixels": round(resolfit, 3),
            "shape_exponent_n": round(float(a[3]), 4),
            "fit_rms": round(rms, 4),
            "center_offset_a0": round(float(a[0]), 3),
            "amplitude_a1": round(float(a[1]), 3),
            "width_a2": round(float(a[2]), 3),
            "peak_pixel_indm": indm,
            "window_lo": lo,
            "window_hi": hi,
        }

        # Full spectrum
        pixel_cols = {f"Pixel_{i}": round(float(signal_mean[i]), 4) for i in range(npix)}
        row = {**summary, **pixel_cols}
        df = pd.DataFrame([row])

        self.output_dir.mkdir(parents=True, exist_ok=True)
        csv_path = self.output_dir / f"CheckSpectrometer_{sn}_{timestamp}.csv"
        df.to_csv(csv_path, index=False)
        LOGGER.info("Check spectrometer CSV saved: %s", csv_path)
        return csv_path
