from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import time
import traceback
import threading
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    import pandas as pd
except Exception:  # pragma: no cover - runtime dependency on target machine
    pd = None

try:
    import serial
except Exception:  # pragma: no cover - runtime dependency on target machine
    serial = None

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

try:
    from characterization_analysis import (
        AnalysisArtifact,
        CharacterizationConfig,
        CharacterizationResult,
        perform_characterization,
    )
except Exception:  # pragma: no cover
    AnalysisArtifact = None
    CharacterizationConfig = None
    CharacterizationResult = None
    perform_characterization = None

try:
    from stage.stage_controller import StageController
except ImportError:
    StageController = None


LOGGER = logging.getLogger(__name__)

KNOWN_HG_AR_NM = np.array(
    [289.36, 296.73, 302.15, 313.16, 334.19, 365.01, 404.66, 407.78, 435.84, 507.30, 546.08],
    dtype=float,
)
IB_REGION_HALF = 20


def get_resource_path(relative_path: str) -> str:
    """Return a bundled-resource path that also works for frozen apps."""
    base_dir = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    return str((base_dir / relative_path).resolve())


def get_writable_path(relative_path: str) -> str:
    """Return a writable path near the app, falling back to the user home."""
    preferred_base = Path(__file__).resolve().parent
    try:
        preferred_base.mkdir(parents=True, exist_ok=True)
        test_path = preferred_base / ".write_test"
        test_path.write_text("ok", encoding="utf-8")
        test_path.unlink(missing_ok=True)
        return str((preferred_base / relative_path).resolve())
    except Exception:
        fallback_base = Path.home() / ".head-scilab"
        fallback_base.mkdir(parents=True, exist_ok=True)
        return str((fallback_base / relative_path).resolve())


def _clean_text(value: object) -> str:
    if isinstance(value, bytes):
        text = value.decode("utf-8", errors="ignore")
    else:
        text = str(value)
    text = text.strip()
    if text.startswith("b'") and text.endswith("'"):
        text = text[2:-1]
    if text.startswith('b"') and text.endswith('"'):
        text = text[2:-1]
    return text.strip()


def _find_peaks(y: Sequence[float], height: Optional[float] = None, distance: int = 1):
    try:
        from scipy.signal import find_peaks  # type: ignore

        return find_peaks(np.asarray(y, dtype=float), height=height, distance=distance)
    except Exception:
        data = np.asarray(y, dtype=float)
        peaks: List[int] = []
        last_idx = -distance
        threshold = -np.inf if height is None else float(height)
        for idx in range(1, len(data) - 1):
            if idx - last_idx < distance:
                continue
            if data[idx] >= threshold and data[idx] > data[idx - 1] and data[idx] >= data[idx + 1]:
                peaks.append(idx)
                last_idx = idx
        return np.asarray(peaks, dtype=int), {"peak_heights": data[peaks] if peaks else np.array([])}


def best_ordered_linear_match(peaks_pix: Sequence[float], candidate_wls: Sequence[float], min_points: int = 5):
    peaks = np.asarray(peaks_pix, dtype=float)
    lines = np.asarray(candidate_wls, dtype=float)
    p_count, l_count = len(peaks), len(lines)
    best = None

    def score(pix_sel: np.ndarray, wl_sel: np.ndarray):
        design = np.vstack([pix_sel, np.ones_like(pix_sel)]).T
        a, b = np.linalg.lstsq(design, wl_sel, rcond=None)[0]
        pred = a * pix_sel + b
        rmse = np.sqrt(np.mean((wl_sel - pred) ** 2))
        return rmse, a, b

    if p_count >= l_count:
        for start in range(p_count - l_count + 1):
            pix_sel = peaks[start : start + l_count]
            wl_sel = lines.copy()
            rmse, a, b = score(pix_sel, wl_sel)
            if best is None or rmse < best[0]:
                best = (rmse, a, b, pix_sel.copy(), wl_sel.copy())
    else:
        for start in range(l_count - p_count + 1):
            pix_sel = peaks.copy()
            wl_sel = lines[start : start + p_count]
            rmse, a, b = score(pix_sel, wl_sel)
            if best is None or rmse < best[0]:
                best = (rmse, a, b, pix_sel.copy(), wl_sel.copy())

    return best if best and len(best[3]) >= min_points else None


def compute_width_at_percent_max(x: Sequence[float], y: Sequence[float], percent: float = 0.5) -> Tuple[float, float, float]:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_arr.size == 0 or y_arr.size == 0 or x_arr.size != y_arr.size:
        return float("nan"), float("nan"), float("nan")

    ymax = float(np.nanmax(y_arr))
    if not np.isfinite(ymax) or ymax <= 0:
        return float("nan"), float("nan"), float("nan")

    level = ymax * float(percent)
    peak_idx = int(np.nanargmax(y_arr))

    left_idx = peak_idx
    while left_idx > 0 and y_arr[left_idx] >= level:
        left_idx -= 1

    right_idx = peak_idx
    while right_idx < len(y_arr) - 1 and y_arr[right_idx] >= level:
        right_idx += 1

    def _interp(idx0: int, idx1: int) -> float:
        x0, y0 = x_arr[idx0], y_arr[idx0]
        x1, y1 = x_arr[idx1], y_arr[idx1]
        if y1 == y0:
            return float(x0)
        return float(x0 + (level - y0) * (x1 - x0) / (y1 - y0))

    x_left = _interp(max(0, left_idx), min(len(x_arr) - 1, left_idx + 1))
    x_right = _interp(max(0, right_idx - 1), min(len(x_arr) - 1, right_idx))
    return x_right - x_left, x_left, x_right


def compute_fwhm(x: Sequence[float], y: Sequence[float]) -> float:
    return compute_width_at_percent_max(x, y, percent=0.5)[0]


def normalized_lsf_from_df(df, wavelength: str, sat_thresh: float = 65400.0):
    if df is None or "Wavelength" not in df.columns:
        return None

    sig_rows = df[df["Wavelength"].astype(str) == str(wavelength)]
    dark_rows = df[df["Wavelength"].astype(str) == f"{wavelength}_dark"]
    if sig_rows.empty:
        return None

    pixel_cols = [col for col in df.columns if str(col).startswith("Pixel_")]
    if not pixel_cols:
        return None

    signal = sig_rows.iloc[-1][pixel_cols].to_numpy(dtype=float)
    dark = np.zeros_like(signal)
    if not dark_rows.empty:
        dark = dark_rows.iloc[-1][pixel_cols].to_numpy(dtype=float)

    if np.nanmax(signal) >= sat_thresh:
        return None

    corrected = signal - dark
    corrected -= np.nanmin(corrected)
    max_val = float(np.nanmax(corrected)) if corrected.size else 0.0
    if not np.isfinite(max_val) or max_val <= 0:
        return None
    return corrected / max_val


def stray_light_metrics(lsf: Sequence[float], peak_pix: int, ib_half: int = 20) -> Dict[str, float]:
    arr = np.asarray(lsf, dtype=float)
    total = float(np.nansum(arr))
    if total <= 0:
        return {"in_band_fraction": np.nan, "out_of_band_fraction": np.nan, "peak_to_total": np.nan}

    lo = max(0, int(peak_pix) - int(ib_half))
    hi = min(len(arr), int(peak_pix) + int(ib_half) + 1)
    in_band = float(np.nansum(arr[lo:hi]))
    out_band = total - in_band
    peak = float(np.nanmax(arr)) if arr.size else np.nan
    return {
        "in_band_fraction": in_band / total,
        "out_of_band_fraction": out_band / total,
        "peak_to_total": peak / total if np.isfinite(peak) else np.nan,
    }


@dataclass
class HardwareState:
    dll_path: str = ""
    spectrometer_type: str = "Auto"
    com_ports: Dict[str, str] = field(default_factory=dict)
    laser_power: Dict[str, float] = field(default_factory=dict)


class MeasurementData:
    def __init__(self, npix: int = 2048, serial_number: str = "Unknown"):
        self.npix = int(npix)
        self.serial_number = serial_number
        self.rows: List[List[float]] = []

    def clear(self) -> None:
        self.rows.clear()

    def to_dataframe(self):
        if pd is None:
            raise RuntimeError("pandas is required to save measurement data.")

        columns = ["Timestamp", "Wavelength", "IntegrationTime", "NumCycles"] + [
            f"Pixel_{idx}" for idx in range(int(self.npix))
        ]

        normalized_rows = []
        for row in self.rows:
            row_values = list(row)
            if len(row_values) < len(columns):
                row_values.extend([np.nan] * (len(columns) - len(row_values)))
            elif len(row_values) > len(columns):
                row_values = row_values[: len(columns)]
            normalized_rows.append(row_values)

        return pd.DataFrame(normalized_rows, columns=columns)

    def last_vectors_for(self, tag: str):
        signal = None
        dark = None
        dark_tag = f"{tag}_dark"

        for row in reversed(self.rows):
            wavelength = str(row[1])
            values = np.asarray(row[4:], dtype=float)
            if signal is None and wavelength == str(tag):
                signal = values
            elif dark is None and wavelength == dark_tag:
                dark = values
            if signal is not None and dark is not None:
                break

        return signal, dark


class SerialDevice:
    def __init__(
        self,
        name: str,
        port: str = "",
        baudrate: int = 9600,
        timeout: float = 1.0,
        line_ending: str = "\r\n",
        **serial_kwargs,
    ):
        self.name = name
        self.port = port
        self.baudrate = baudrate
        self.timeout = timeout
        self.line_ending = line_ending
        self.serial_kwargs = serial_kwargs
        self._ser = None

    def configure(self, port: Optional[str] = None, baudrate: Optional[int] = None, **serial_kwargs) -> None:
        if port is not None:
            self.port = port
        if baudrate is not None:
            self.baudrate = baudrate
        if serial_kwargs:
            self.serial_kwargs.update(serial_kwargs)

    def open(self) -> bool:
        if not self.port or serial is None:
            return False
        if self._ser is not None and getattr(self._ser, "is_open", False):
            if getattr(self._ser, "port", None) == self.port:
                return True
            self.close()

        try:
            self._ser = serial.Serial(
                self.port,
                baudrate=self.baudrate,
                timeout=self.timeout,
                write_timeout=self.timeout,
                **self.serial_kwargs,
            )
            return True
        except Exception:
            LOGGER.exception("Unable to open serial port %s for %s", self.port, self.name)
            self._ser = None
            return False

    def close(self) -> None:
        if self._ser is not None:
            try:
                self._ser.close()
            except Exception:
                pass
            self._ser = None

    def reset_buffers(self) -> None:
        if self._ser is None:
            return
        for method_name in ("reset_input_buffer", "reset_output_buffer"):
            try:
                getattr(self._ser, method_name)()
            except Exception:
                pass

    def write_text(self, text: str) -> None:
        if self._ser is None:
            raise RuntimeError(f"{self.name} serial port is not open.")
        payload = (text + self.line_ending).encode("utf-8")
        self._ser.write(payload)
        self._ser.flush()

    def read_all_text(self) -> str:
        if self._ser is None:
            return ""
        try:
            data = self._ser.read_all()
        except Exception:
            return ""
        return data.decode("utf-8", errors="ignore").strip()


class LaserController:
    OBIS_MAP = {"405": 5, "445": 4, "488": 3, "640": 2, "685": 6}

    def __init__(self, com_ports: Optional[Dict[str, str]] = None):
        cube_kwargs = {}
        if serial is not None:
            cube_kwargs = {
                "parity": serial.PARITY_NONE,
                "stopbits": serial.STOPBITS_ONE,
                "bytesize": serial.EIGHTBITS,
            }

        self.obis = SerialDevice("OBIS", baudrate=9600, timeout=1.0, line_ending="\r\n")
        self.cube = SerialDevice("CUBE", baudrate=19200, timeout=1.0, line_ending="\r", **cube_kwargs)
        self.relay = SerialDevice("RELAY", baudrate=9600, timeout=1.0, line_ending="\r")
        self.configure_ports(com_ports or {})

    def configure_ports(self, com_ports: Dict[str, str]) -> None:
        self.obis.configure(port=com_ports.get("OBIS", ""))
        self.cube.configure(port=com_ports.get("CUBE", ""))
        self.relay.configure(port=com_ports.get("RELAY", ""))

    def open_all(self) -> None:
        for device in (self.obis, self.cube, self.relay):
            if device.port:
                device.open()

    def ensure_open_for_tag(self, tag: str) -> None:
        if tag in self.OBIS_MAP:
            if not self.obis.open():
                raise RuntimeError(f"Unable to open OBIS port '{self.obis.port}'.")
        elif tag == "377":
            if not self.cube.open():
                raise RuntimeError(f"Unable to open CUBE port '{self.cube.port}'.")
        elif tag in {"517", "532", "Hg_Ar"}:
            if not self.relay.open():
                raise RuntimeError(f"Unable to open RELAY port '{self.relay.port}'.")

    def obis_cmd(self, command: str) -> str:
        if not self.obis.open():
            raise RuntimeError(f"Unable to open OBIS port '{self.obis.port}'.")
        self.obis.reset_buffers()
        self.obis.write_text(command)
        time.sleep(0.2)
        return self.obis.read_all_text()

    def obis_on(self, channel: int) -> None:
        self.obis_cmd(f"SOUR{channel}:AM:STAT ON")

    def obis_off(self, channel: int) -> None:
        self.obis_cmd(f"SOUR{channel}:AM:STAT OFF")

    def obis_set_power(self, channel: int, power_watts: float) -> None:
        self.obis_cmd(f"SOUR{channel}:POW:LEV:IMM:AMPL {float(power_watts):.4f}")

    def cube_cmd(self, command: str) -> str:
        if not self.cube.open():
            raise RuntimeError(f"Unable to open CUBE port '{self.cube.port}'.")
        last_response = ""
        for _ in range(3):
            self.cube.reset_buffers()
            self.cube.write_text(command)
            time.sleep(0.4)
            last_response = self.cube.read_all_text()
            if last_response:
                break
        return last_response

    def cube_on(self, power_mw: float = 12.0) -> None:
        self.cube_cmd("EXT=1")
        time.sleep(1.0)
        self.cube_cmd("CW=1")
        time.sleep(1.0)
        self.cube_cmd(f"P={float(power_mw):.3f}")
        time.sleep(1.0)
        self.cube_cmd("L=1")
        time.sleep(3.0)

    def cube_off(self) -> None:
        self.cube_cmd("L=0")

    def relay_cmd(self, command: str) -> str:
        if not self.relay.open():
            raise RuntimeError(f"Unable to open RELAY port '{self.relay.port}'.")
        self.relay.write_text(command)
        time.sleep(0.05)
        return self.relay.read_all_text()

    def relay_on(self, relay_number: int) -> None:
        self.relay_cmd(f"R{int(relay_number)}S")

    def relay_off(self, relay_number: int) -> None:
        self.relay_cmd(f"R{int(relay_number)}R")

    def all_off(self) -> None:
        for channel in self.OBIS_MAP.values():
            try:
                self.obis_off(channel)
            except Exception:
                pass
        for relay_number in (1, 2, 3, 4):
            try:
                self.relay_off(relay_number)
            except Exception:
                pass
        try:
            self.cube_off()
        except Exception:
            pass


class FilterWheelController:
    def __init__(self, port: str = ""):
        self.device = SerialDevice("HEADSENSOR", port=port, baudrate=4800, timeout=1.0, line_ending="\r")
        self.serial_status = {"hst": ["Head sensor idle."]}

    def configure_port(self, port: str) -> None:
        self.device.configure(port=port)

    def _record_status(self, status: str) -> None:
        messages = self.serial_status.setdefault("hst", [])
        messages.append(status)
        if len(messages) > 50:
            del messages[:-50]

    def open(self) -> bool:
        ok = self.device.open()
        self._record_status("opened" if ok else "open failed")
        return ok

    def close(self) -> None:
        self.device.close()
        self._record_status("closed")

    def send_raw_command(self, command: str, pause_s: float = 0.2) -> Tuple[bool, str]:
        if not self.device.open():
            status = f"Unable to open head sensor port '{self.device.port}'."
            self._record_status(status)
            return False, status

        try:
            self.device.reset_buffers()
            self.device.write_text(command)
            time.sleep(pause_s)
            response = self.device.read_all_text()
            status = response or "timeout"
            self._record_status(status)
            return True, status
        except Exception as exc:
            status = str(exc)
            self._record_status(status)
            return False, status

    def query_device_id(self) -> Tuple[bool, str]:
        success, response = self.send_raw_command("?", pause_s=0.3)
        cleaned = _clean_text(response)
        if cleaned.startswith("Pan"):
            return True, cleaned
        return success and bool(cleaned), cleaned or response

    def set_filterwheel(self, fw_num: int, position: int) -> bool:
        success, response = self.send_raw_command(f"F{int(fw_num)}{int(position)}", pause_s=0.3)
        return success and "error" not in response.lower()

    def reset_filterwheel(self, fw_num: int) -> bool:
        success, response = self.send_raw_command(f"F{int(fw_num)}r", pause_s=0.5)
        return success and "error" not in response.lower()

    def test_filterwheel(self, fw_num: int) -> bool:
        success, response = self.send_raw_command(f"F{int(fw_num)}m", pause_s=0.3)
        return success and "error" not in response.lower()


class SpectroApp(tk.Tk):
    DEFAULT_SPECTROMETER_TYPE = "Auto"
    DEFAULT_ALL_LASERS = ["377", "405", "445", "488", "517", "532", "640", "685", "Hg_Ar"]
    DEFAULT_COM_PORTS = {"OBIS": "COM10", "CUBE": "COM1", "RELAY": "COM11", "HEADSENSOR": "", "STAGE": ""}
    DEFAULT_LASER_POWERS = {
        "377": 12.0,
        "405": 0.05,
        "445": 0.03,
        "488": 0.03,
        "517": 30.0,
        "532": 30.0,
        "640": 0.05,
        "685": 0.03,
        "Hg_Ar": 0.0,
    }
    DEFAULT_START_IT = {"532": 20.0, "517": 80.0, "Hg_Ar": 20.0, "default": 2.4}

    N_SIG = 50
    N_DARK = 50
    N_SIG_640 = 10
    N_DARK_640 = 10

    TARGET_LOW = 60000
    TARGET_HIGH = 65000
    TARGET_MID = 62500

    IT_MIN = 0.2
    IT_MAX = 3000.0
    IT_STEP_UP = 0.3
    IT_STEP_DOWN = 0.1
    MAX_IT_ADJUST_ITERS = 1000
    SAT_THRESH = 65400

    SETTINGS_FILE = get_writable_path("spectro_gui_settings.json")
    RESULTS_ROOT = Path(get_writable_path("results"))

    def _configure_initial_geometry(self) -> None:
        screen_width = max(800, self.winfo_screenwidth())
        screen_height = max(600, self.winfo_screenheight())

        target_width = min(1500, int(screen_width * 0.96), max(820, screen_width - 80))
        target_height = min(950, int(screen_height * 0.92), max(620, screen_height - 80))

        min_width = min(target_width, 960 if screen_width >= 1100 else max(720, int(screen_width * 0.78)))
        min_height = min(target_height, 680 if screen_height >= 760 else max(540, int(screen_height * 0.8)))

        self.minsize(min_width, min_height)

        offset_x = max(0, (screen_width - target_width) // 2)
        offset_y = max(0, (screen_height - target_height) // 2)
        self.geometry(f"{target_width}x{target_height}+{offset_x}+{offset_y}")

    def __init__(self):
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
        super().__init__()

        self.title("SciGlob Spectrometer Characterization System")
        self._configure_initial_geometry()

        try:
            self.iconbitmap(get_resource_path("sciglob_symbol.ico"))
        except Exception:
            pass

        self.npix = 2048
        self.sn = "Unknown"
        self.spec = None
        self.spec_backend = None

        self.hw = HardwareState(
            dll_path="",
            spectrometer_type=self.DEFAULT_SPECTROMETER_TYPE,
            com_ports=dict(self.DEFAULT_COM_PORTS),
            laser_power=dict(self.DEFAULT_LASER_POWERS),
        )
        self.data = MeasurementData(npix=self.npix, serial_number=self.sn)
        self.available_lasers = list(self.DEFAULT_ALL_LASERS)
        self.laser_configs = self._build_default_laser_configs()

        self.lasers = LaserController(self.hw.com_ports)
        self.filterwheel = FilterWheelController(self.hw.com_ports["HEADSENSOR"])

        self.stage = StageController() if StageController else None
        self.stage_config_path = ""

        self.live_running = threading.Event()
        self.measure_running = threading.Event()
        self._pending_it = None
        self._it_updating = False
        self.it_history: List[Tuple[float, float]] = []

        self.analysis_artifacts: List = []
        self.analysis_summary_lines: List[str] = []
        self.analysis_measurement_tabs: Dict[str, object] = {}
        self.analysis_measurement_counter = 0
        self.reference_csv_paths: List[str] = []
        self._latest_results_dir: Optional[Path] = None
        self._latest_csv_path: Optional[Path] = None
        self._latest_results_timestamp: Optional[str] = None
        self._analysis_images: List[object] = []

        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True)

        self.setup_tab = ttk.Frame(self.nb)
        self.live_tab = ttk.Frame(self.nb)
        self.measure_tab = ttk.Frame(self.nb)
        self.analysis_tab = ttk.Frame(self.nb)
        self.eeprom_tab = ttk.Frame(self.nb)

        self.nb.add(self.setup_tab, text="Setup")
        self.nb.add(self.live_tab, text="Live View")
        self.nb.add(self.measure_tab, text="Measurements")
        self.nb.add(self.analysis_tab, text="Analysis")
        self.nb.add(self.eeprom_tab, text="EEPROM")

        from tabs import analysis_tab, eeprom_tab, live_view_tab, measurements_tab, setup_tab

        setup_tab.build(self)
        live_view_tab.build(self)
        measurements_tab.build(self)
        analysis_tab.build(self)
        eeprom_tab.build(self)

        if hasattr(self, "load_settings_into_ui"):
            try:
                self.load_settings_into_ui()
            except Exception:
                LOGGER.exception("Failed to load saved settings")

        self.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _build_default_laser_configs(self) -> Dict[str, Dict[str, float]]:
        configs = {}
        for laser in self.DEFAULT_ALL_LASERS:
            if laser == "377":
                laser_type = "CUBE"
            elif laser in {"517", "532", "Hg_Ar"}:
                laser_type = "RELAY"
            else:
                laser_type = "OBIS"
            configs[laser] = {"type": laser_type, "power": self.DEFAULT_LASER_POWERS.get(laser, 0.01)}
        return configs

    def _post_error(self, title: str, ex: Exception) -> None:
        tb = "".join(traceback.format_exception(type(ex), ex, ex.__traceback__))
        LOGGER.error("[%s] %s\n%s", title, ex, tb)
        try:
            self.after(0, lambda: messagebox.showerror(title, str(ex)))
        except Exception:
            pass

    def _on_closing(self) -> None:
        handler = getattr(self, "on_close", None)
        if callable(handler):
            try:
                handler()
                return
            except Exception:
                LOGGER.exception("Error while closing application")
        self.destroy()

    def _live_reset_view(self) -> None:
        try:
            self.live_limits_locked = False
            if hasattr(self, "live_ax"):
                self.live_ax.relim()
                self.live_ax.autoscale_view()
            if hasattr(self, "live_canvas"):
                self.live_canvas.draw_idle()
        except Exception as exc:
            self._post_error("Reset Zoom", exc)

    def rebuild_laser_ui(self) -> None:
        self.available_lasers = list(self.laser_configs.keys())
        if "Hg_Ar" not in self.available_lasers:
            self.available_lasers.append("Hg_Ar")

    def update_target_peak(self, value) -> None:
        try:
            target = int(value)
        except Exception:
            return
        self.TARGET_MID = target
        self.TARGET_LOW = max(0, target - 2500)
        self.TARGET_HIGH = target + 2500
        if hasattr(self, "target_band_label"):
            self.target_band_label.config(text=f"Target window: {self.TARGET_LOW}-{self.TARGET_HIGH}")

    def _clear_analysis_window(self) -> None:
        """Reset current-run analysis state. Previous run tabs are preserved."""
        self.analysis_artifacts = []
        self.analysis_summary_lines = []
        self._analysis_images = []
        self._latest_results_dir = None
        self._latest_csv_path = None
        self._latest_results_timestamp = None

    def _ensure_results_dir(self) -> Path:
        self.RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
        if self._latest_results_dir is None:
            stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            serial_number = _clean_text(self.data.serial_number or self.sn or "Unknown")
            serial_number = serial_number.replace(" ", "_")
            self._latest_results_dir = self.RESULTS_ROOT / f"{serial_number}_{stamp}"
            self._latest_results_dir.mkdir(parents=True, exist_ok=True)
        return self._latest_results_dir

    def save_measurement_data(self) -> Optional[str]:
        if not self.data.rows:
            return None
        if pd is None:
            raise RuntimeError("pandas is required to save measurement data.")

        output_dir = self._ensure_results_dir()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        serial_number = _clean_text(self.data.serial_number or self.sn or "Unknown")
        path = output_dir / f"Measurements_{serial_number}_{timestamp}.csv"
        df = self.data.to_dataframe()
        df.to_csv(path, index=False)
        self._latest_csv_path = path
        return str(path)

    def _figure_to_artifact(self, title: str, figure: Figure, output_path: Path) -> Dict[str, object]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=150, bbox_inches="tight")
        return {"title": title, "figure": figure, "path": str(output_path)}

    def run_analysis_and_save_plots(self, csv_path: Optional[str] = None):
        """Run comprehensive analysis using perform_characterization and save all plots."""
        if pd is None:
            raise RuntimeError("pandas is required to run analysis.")
        if perform_characterization is None:
            raise RuntimeError("characterization_analysis module is required for analysis.")

        if csv_path:
            df = pd.read_csv(csv_path)
            results_dir = Path(csv_path).resolve().parent
            self._latest_csv_path = Path(csv_path).resolve()
            self._latest_results_dir = results_dir
        else:
            if not self.data.rows:
                return []
            df = self.data.to_dataframe()
            results_dir = self._ensure_results_dir()
            if self._latest_csv_path is None:
                saved_path = self.save_measurement_data()
                if saved_path:
                    self._latest_csv_path = Path(saved_path)

        plots_dir = results_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        if "Wavelength" in df.columns:
            df["Wavelength"] = df["Wavelength"].astype(str)

        sn = _clean_text(self.data.serial_number or self.sn or "Unknown")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        reference_csv_paths = getattr(self, "reference_csv_paths", [])

        LOGGER.info("Starting comprehensive analysis via perform_characterization...")
        try:
            result = perform_characterization(
                df, sn, str(plots_dir), timestamp, None, reference_csv_paths or None
            )
        except Exception as exc:
            LOGGER.exception("perform_characterization raised an exception")
            self.analysis_artifacts = []
            self.analysis_summary_lines = [f"Analysis error: {exc}"]
            self._latest_results_timestamp = timestamp
            return []

        self.analysis_artifacts = result.artifacts if result else []
        self.analysis_summary_lines = result.summary_lines if result else []
        self._latest_results_timestamp = timestamp

        if hasattr(self, "export_plots_btn"):
            self.export_plots_btn.state(["!disabled"] if self.analysis_artifacts else ["disabled"])
        if hasattr(self, "open_folder_btn"):
            self.open_folder_btn.state(["!disabled"] if self.analysis_artifacts else ["disabled"])

        paths = [getattr(art, "path", "") for art in self.analysis_artifacts]
        LOGGER.info("Analysis complete. %d plots generated to: %s", len(paths), plots_dir)
        return paths

    def refresh_analysis_view(self):
        """Add a new run tab to the analysis notebook with all plots from the latest analysis."""
        if not self.analysis_artifacts:
            if self._latest_csv_path and Path(self._latest_csv_path).is_file():
                self.run_analysis_and_save_plots(str(self._latest_csv_path))
            elif self.data.rows:
                csv_path = self.save_measurement_data()
                if csv_path:
                    self.run_analysis_and_save_plots(csv_path)
            else:
                messagebox.showwarning("Analysis", "No measurement data available.")
                return

        if not self.analysis_artifacts:
            messagebox.showwarning("Analysis", "No analysis plots were generated.")
            return

        notebook = getattr(self, "analysis_measurements_notebook", None)
        if notebook is None:
            return

        self.analysis_measurement_counter += 1
        timestamp = getattr(self, "_latest_results_timestamp", None) or datetime.now().strftime("%Y%m%d_%H%M%S")
        tab_name = f"Run #{self.analysis_measurement_counter} ({timestamp})"

        # Remove welcome tab if it exists
        if hasattr(self, "analysis_welcome_tab"):
            try:
                notebook.forget(self.analysis_welcome_tab)
            except Exception:
                pass

        # Create new tab for this measurement run
        measurement_tab = ttk.Frame(notebook)
        notebook.add(measurement_tab, text=tab_name)

        # Scrollable canvas for plots
        canvas_container = tk.Canvas(measurement_tab, highlightthickness=0)
        scrollbar_y = ttk.Scrollbar(measurement_tab, orient="vertical", command=canvas_container.yview)
        scrollbar_x = ttk.Scrollbar(measurement_tab, orient="horizontal", command=canvas_container.xview)
        scrollable_frame = ttk.Frame(canvas_container)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas_container.configure(scrollregion=canvas_container.bbox("all")),
        )
        canvas_container.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas_container.configure(yscrollcommand=scrollbar_y.set, xscrollcommand=scrollbar_x.set)

        def _on_mousewheel(event):
            canvas_container.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas_container.bind("<Enter>", lambda e: canvas_container.bind_all("<MouseWheel>", _on_mousewheel))
        canvas_container.bind("<Leave>", lambda e: canvas_container.unbind_all("<MouseWheel>"))

        canvas_container.grid(row=0, column=0, sticky="nsew")
        scrollbar_y.grid(row=0, column=1, sticky="ns")
        scrollbar_x.grid(row=1, column=0, sticky="ew")
        measurement_tab.grid_rowconfigure(0, weight=1)
        measurement_tab.grid_columnconfigure(0, weight=1)

        # Arrange plots in a 2-column grid
        num_plots = len(self.analysis_artifacts)
        cols = 2 if num_plots > 1 else 1

        for idx, artifact in enumerate(self.analysis_artifacts):
            row = idx // cols
            col = idx % cols

            art_name = getattr(artifact, "name", f"Plot {idx + 1}")
            art_fig = getattr(artifact, "figure", None)
            if art_fig is None:
                continue

            plot_frame = ttk.LabelFrame(scrollable_frame, text=art_name, padding=12)
            plot_frame.grid(row=row, column=col, padx=10, pady=10, sticky="nsew")

            try:
                art_fig.set_size_inches(8, 6)
                art_fig.set_dpi(100)
                art_fig.tight_layout(pad=1.5)
                for ax in art_fig.get_axes():
                    ax.grid(True, alpha=0.2, linestyle="--", linewidth=0.5)
                art_fig.patch.set_facecolor("white")
            except Exception:
                pass

            try:
                fig_canvas = FigureCanvasTkAgg(art_fig, master=plot_frame)
                fig_canvas.draw()
                canvas_widget = fig_canvas.get_tk_widget()
                canvas_widget.configure(highlightthickness=0)
                canvas_widget.pack(fill="both", expand=True, padx=5, pady=5)
            except Exception:
                ttk.Label(plot_frame, text="Plot rendering failed.",
                          foreground="red").pack(padx=5, pady=5)
                continue

            try:
                toolbar_frame = ttk.Frame(plot_frame)
                toolbar_frame.pack(fill="x", padx=5, pady=(0, 5))
                toolbar = NavigationToolbar2Tk(fig_canvas, toolbar_frame)
                toolbar.update()
            except Exception:
                pass

        for c in range(cols):
            scrollable_frame.grid_columnconfigure(c, weight=1, uniform="plots")

        # Summary section at the bottom
        summary_row = (num_plots + cols - 1) // cols
        summary_frame = ttk.LabelFrame(scrollable_frame, text="Analysis Summary", padding=12)
        summary_frame.grid(row=summary_row, column=0, columnspan=cols, padx=10, pady=10, sticky="ew")

        summary_text_widget = tk.Text(
            summary_frame, height=10, wrap="word",
            font=("TkDefaultFont", 9), relief="flat",
            bg="#f8f9fa", fg="#212529", padx=10, pady=10,
        )
        summary_text = "\n".join(self.analysis_summary_lines) if self.analysis_summary_lines else "No summary available."
        summary_text_widget.insert("1.0", summary_text)
        summary_text_widget.configure(state="disabled")
        summary_text_widget.pack(fill="both", expand=True)

        # Store tab reference
        self.analysis_measurement_tabs[tab_name] = {
            "tab": measurement_tab,
            "timestamp": timestamp,
            "csv_path": str(self._latest_csv_path) if self._latest_csv_path else None,
        }

        if hasattr(self, "export_plots_btn"):
            self.export_plots_btn.state(["!disabled"])
        if hasattr(self, "open_folder_btn"):
            self.open_folder_btn.state(["!disabled"])

        # Switch to the new tab and the Analysis page
        notebook.select(measurement_tab)
        try:
            self.nb.select(self.analysis_tab)
        except Exception:
            pass

    def export_analysis_plots(self):
        if not self.analysis_artifacts:
            messagebox.showwarning("Export Plots", "No analysis plots available.")
            return
        destination = filedialog.askdirectory(title="Select folder to export analysis plots")
        if not destination:
            return
        dest_dir = Path(destination)
        copied = 0
        for artifact in self.analysis_artifacts:
            src = Path(str(getattr(artifact, "path", "")))
            if not src.is_file():
                continue
            shutil.copy2(src, dest_dir / src.name)
            copied += 1
        messagebox.showinfo("Export Plots", f"Exported {copied} plot(s) to:\n{dest_dir}")

    def export_analysis_summary(self):
        if not self.analysis_summary_lines:
            messagebox.showwarning("Export Summary", "No analysis summary available.")
            return
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = filedialog.asksaveasfilename(
            title="Save Analysis Summary",
            defaultextension=".txt",
            initialfile=f"analysis_summary_{timestamp}.txt",
            filetypes=[("Text", "*.txt"), ("All Files", "*.*")],
        )
        if not path:
            return
        Path(path).write_text("\n".join(self.analysis_summary_lines) + "\n", encoding="utf-8")
        messagebox.showinfo("Export Summary", f"Saved summary to:\n{path}")

    def open_results_folder(self):
        if self._latest_results_dir is None or not self._latest_results_dir.exists():
            messagebox.showwarning("Results Folder", "No results folder available yet.")
            return
        folder = self._latest_results_dir.resolve()
        if not folder.is_dir():
            messagebox.showwarning("Results Folder", "Path is not a valid directory.")
            return
        results_root = self.RESULTS_ROOT.resolve()
        if not str(folder).startswith(str(results_root)):
            messagebox.showerror("Results Folder", "Path is outside the results directory.")
            return
        try:
            if sys.platform.startswith("win"):
                os.startfile(str(folder))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.Popen(["open", "--", str(folder)])  # noqa: S603
            else:
                subprocess.Popen(["xdg-open", "--", str(folder)])  # noqa: S603
        except Exception as exc:
            self._post_error("Open Results Folder", exc)

    def _update_auto_it_plot(self, tag, spectrum, it_ms, peak):
        """Update measurement plot during Auto-IT adjustment (called from UI thread)."""
        try:
            xs = np.arange(len(spectrum))
            self.meas_sig_line.set_data(xs, spectrum)
            self.meas_ax.set_xlim(0, max(10, len(spectrum) - 1))
            self.meas_ax.set_ylim(0, 65000)
            self.meas_ax.set_title(
                f"Auto-IT: {tag} nm | IT={it_ms:.2f} ms | Peak={peak:.0f}",
                fontsize=13, fontweight="bold", pad=8,
            )

            steps = list(getattr(self, "it_history", []))
            if steps:
                st = np.arange(len(steps))
                peaks = [p for (_, p) in steps]
                its = [it for (it, _) in steps]
                self.inset_peak_line.set_data(st, peaks)
                self.inset_it_line.set_data(st, its)
                self.meas_inset.set_xlim(-0.5, max(0.5, len(st) - 0.5))
                self.meas_inset.relim()
                self.meas_inset.autoscale_view()
                self.meas_inset2.relim()
                self.meas_inset2.autoscale_view()

            self.meas_canvas.draw_idle()
        except Exception:
            LOGGER.debug("Auto-IT plot update skipped (UI not ready)")

    def _finalize_measurement_run(self):
        return None
