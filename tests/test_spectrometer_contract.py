from __future__ import annotations

import unittest

import numpy as np

from domain.spectrometer import assert_spectrometer_backend, validate_spectrometer_backend


class ValidBackend:
    sn = "VALID-001"
    npix_active = 4
    rcm = np.zeros(4, dtype=float)
    abort_on_saturation = False

    def connect(self):
        return "OK"

    def disconnect(self, *args, **kwargs):
        return "OK"

    def set_it(self, it_ms: float):
        return "OK"

    def measure(self, ncy: int = 1):
        return "OK"

    def wait_for_measurement(self):
        return "OK"


class InvalidBackend:
    sn = "BROKEN-001"


class SpectrometerContractTests(unittest.TestCase):
    def test_accepts_valid_backend_shape(self):
        backend = assert_spectrometer_backend(ValidBackend())
        self.assertEqual(backend.sn, "VALID-001")

    def test_reports_missing_members(self):
        issues = validate_spectrometer_backend(InvalidBackend())
        self.assertTrue(any("npix_active" in issue for issue in issues))
        with self.assertRaises(TypeError):
            assert_spectrometer_backend(InvalidBackend())


if __name__ == "__main__":
    unittest.main()

