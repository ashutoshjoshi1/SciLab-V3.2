from __future__ import annotations

import unittest

import numpy as np

from domain.measurement import MeasurementData
from services.measurement_orchestrator import (
    MeasurementOrchestrator,
    MeasurementOrchestratorCallbacks,
    MeasurementOrchestratorConfig,
)


class FakeLasers:
    OBIS_MAP = {"405": 5, "640": 2}

    def __init__(self):
        self.active_tag = None

    def open_all(self):
        return None

    def all_off(self):
        self.active_tag = None

    def ensure_open_for_tag(self, _tag: str):
        return None

    def obis_set_power(self, _channel: int, _power: float):
        return None

    def obis_on(self, channel: int):
        if channel == 5:
            self.active_tag = "405"
        elif channel == 2:
            self.active_tag = "640"

    def obis_off(self, channel: int):
        if (channel == 5 and self.active_tag == "405") or (channel == 2 and self.active_tag == "640"):
            self.active_tag = None

    def cube_on(self, power_mw: float = 0.0):
        self.active_tag = "377"

    def cube_off(self):
        if self.active_tag == "377":
            self.active_tag = None

    def relay_on(self, relay_number: int):
        self.active_tag = {1: "532", 2: "517", 4: "Hg_Ar"}.get(relay_number)

    def relay_off(self, relay_number: int):
        if self.active_tag == {1: "532", 2: "517", 4: "Hg_Ar"}.get(relay_number):
            self.active_tag = None


class FakeSpectrometer:
    def __init__(self, lasers: FakeLasers):
        self.lasers = lasers
        self.sn = "SN-FAKE"
        self.npix_active = 8
        self.rcm = np.zeros(self.npix_active, dtype=float)
        self.abort_on_saturation = False
        self._integration_time_ms = 1.0

    def connect(self):
        return "OK"

    def disconnect(self, *args, **kwargs):
        return "OK"

    def set_it(self, it_ms: float):
        self._integration_time_ms = float(it_ms)
        return "OK"

    def measure(self, ncy: int = 1):
        return "OK"

    def wait_for_measurement(self):
        if self.lasers.active_tag is None:
            self.rcm = np.full(self.npix_active, 50.0, dtype=float)
        else:
            self.rcm = np.full(self.npix_active, self._integration_time_ms * 100.0, dtype=float)
        return "OK"


class MeasurementOrchestratorTests(unittest.TestCase):
    def test_records_signal_and_dark_measurements_for_standard_run(self):
        lasers = FakeLasers()
        spectrometer = FakeSpectrometer(lasers)
        data = MeasurementData(npix=8, serial_number="SN-FAKE")
        config = MeasurementOrchestratorConfig(
            default_start_it={"default": 10.0, "405": 10.0},
            target_low=900.0,
            target_high=1100.0,
            target_mid=1000.0,
            it_min=0.1,
            it_max=50.0,
            it_step_up=1.0,
            it_step_down=1.0,
            max_it_adjust_iters=5,
            sat_thresh=65000.0,
            n_sig=2,
            n_dark=1,
            n_sig_640=1,
            n_dark_640=1,
        )
        callbacks = MeasurementOrchestratorCallbacks(
            prepare_devices=lambda: None,
            power_lookup=lambda _tag: 0.1,
        )
        orchestrator = MeasurementOrchestrator(
            spectrometer=spectrometer,
            laser_controller=lasers,
            measurement_data=data,
            config=config,
            callbacks=callbacks,
            sleep_fn=lambda _seconds: None,
        )

        result = orchestrator.run(["405"], should_continue=lambda: True)

        self.assertEqual(result.completed_tags, ["405"])
        self.assertEqual(result.rows_written, 2)
        self.assertEqual(len(data.rows), 2)
        self.assertEqual(data.rows[0][1], "405")
        self.assertEqual(data.rows[1][1], "405_dark")


if __name__ == "__main__":
    unittest.main()

