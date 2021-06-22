
import numpy as np
import unittest
import sys
import os
from workflow.scripts.calc_helpers import calc_sample_standard_deviation
from workflow.scripts.calc_helpers import calc_95_confidence_interval
from workflow.scripts.calc_helpers import calc_cleavage_fraction
from workflow.scripts.calc_helpers import calc_fold_change

# set module path for testing
# sys.path.insert(0, "path_in_PYTHONPATH")


class TestSuiteCalcHelpers(unittest.TestCase):
    """Test suite to test the calc helpers functions"""

    def test_calc_sample_standard_deviation(self):
        """Function tests the calc_sample_standard_deviation() function."""
        test_array_1 = np.array(
            [9, 2, 5, 4, 1, 7, 8, 11, 9, 3, 7, 4, 12, 5, 4, 10, 9, 6, 9, 4], dtype=np.float32)
        test_array_2 = np.array(
            [4, 9, 11, 12, 17, 5, 8, 12, 14], dtype=np.float32)
        self.assertAlmostEqual(
            calc_sample_standard_deviation(test_array_1), 3.10305449)
        self.assertAlmostEqual(
            calc_sample_standard_deviation(test_array_2), 4.17665472)

    def test_calc_95_confidence_interval(self):
        """Function tests the calc_95_confidence_interval() function."""
        self.assertEqual(
            type(calc_95_confidence_interval(5, 2)), type((1.08, 8.92)), "Function does not return a tuple")
        self.assertEqual(
            len(calc_95_confidence_interval(5, 2)), len((1.08, 8.92)), "Function does not return a tuple with the correct length")
        self.assertSequenceEqual(
            calc_95_confidence_interval(5, 2), (1.08, 8.92))
        self.assertSequenceEqual(
            calc_95_confidence_interval(0, 3), (-5.88, 5.88))
        self.assertSequenceEqual(
            calc_95_confidence_interval(0, 0), (0.0, 0.0))

    def test_calc_cleavage_fraction(self):
        pass

    def test_calc_fold_change(self):
        pass
