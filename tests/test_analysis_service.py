from __future__ import annotations

import tempfile
import unittest

import pandas as pd

from services.analysis_service import AnalysisService


class AnalysisServiceTests(unittest.TestCase):
    def test_handles_missing_pixel_columns_gracefully(self):
        dataframe = pd.DataFrame(
            {
                "Wavelength": ["405"],
                "IntegrationTime": [1.0],
                "NumCycles": [1],
            }
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            result = AnalysisService().analyze(
                df=dataframe,
                serial_number="SN-TEST",
                output_dir=temp_dir,
                timestamp="20260406_120000",
            )

        self.assertEqual(result.artifacts, [])
        self.assertIn("No pixel columns detected", result.summary_text)


if __name__ == "__main__":
    unittest.main()

