import os
import sys
import unittest
import numpy as np
import workflow.scripts.calc_helpers

# set module path for testing
# sys.path.insert(0, "path_in_PYTHONPATH")


class TestSuiteCalcHelpers(unittest.TestCase):
    """Test suite to test the calc helpers functions"""

    def test_calculate_sample_standard_deviation(self):
        test_array = np.array(
            [9, 2, 5, 4, 1, 7, 8, 11, 9, 3, 7, 4, 12, 5, 4, 10, 9, 6, 9, 4], dtype=np.float32)
        self.assertAlmostEqual(
            calc_helpers.calculate_sample_standard_deviation(test_array), 3.061)
