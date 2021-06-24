import unittest
import regex
from workflow.scripts.yaml_read_helpers import retrieve_dataset_name
from workflow.scripts.yaml_read_helpers import retrieve_minimum_reads
from workflow.scripts.yaml_read_helpers import retrieve_compiled_patterns
from workflow.scripts.yaml_read_helpers import retrieve_prefix_name
from workflow.scripts.yaml_read_helpers import retrieve_compiled_info_patterns


class TestSuiteCalcHelpers(unittest.TestCase):
    """Test suite to test the calc helpers functions"""

    def test_retrieve_dataset_name(self):

        yaml_file_path = "./config/test_config.yaml"

        self.assertEqual(retrieve_dataset_name(yaml_file_path), "T1_D80")

    def test_retrieve_minimum_reads(self):
        yaml_file_path = "./config/test_config.yaml"

        self.assertEqual(retrieve_minimum_reads(yaml_file_path), 30)

    def test_retrieve_compiled_patterns(self):
        """Function tests the retrieve_compiled_patterns() function."""

        yaml_file_path = "./config/test_config.yaml"

        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="prefix", ligand_present=True),
            {regex.compile("(?e)(AGATCTTTTCCGTATATCTCGCCAG){e<=1}"): "A",
             regex.compile("(?e)(AGATGGGAAACAAACAAA){e<=1}"): "W",
             regex.compile("(?e)(AGATACAAAACAAAAC){e<=1}"): "Z"})
        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="suffix", ligand_present=True),
            {regex.compile("(?e)(AAAAAGAAATT){e<=1}"): "X"})
        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="prefix", ligand_present=False),
            {regex.compile("(?e)(CTTTTCCGTATATCTCGCCAG){e<=1}"): "A",
             regex.compile("(?e)(GGGAAACAAACAAA){e<=1}"): "W",
             regex.compile("(?e)(ACAAAACAAAAC){e<=1}"): "Z"})
        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="suffix", ligand_present=False),
            {regex.compile("(?e)(AAAAAGAAAT){e<=1}"): "X"})

    def test_retrieve_prefix_name(self):
        """Function tests the retrieve_prefix_name() function."""

        yaml_file_path = "./config/test_config.yaml"

        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=True, ligand_present=True), "A")
        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=False, ligand_present=True), "W")
        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=True, ligand_present=False), "A")
        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=False, ligand_present=False), "W")

    def test_retrieve_compiled_info_patterns(self):
        """Function tests the retrieve_compiled_info_patterns() function."""

        yaml_file_path = "./config/test_config.yaml"

        self.assertEqual(retrieve_compiled_info_patterns(yaml_file_path), (regex.compile(
            r"(?<=D)([0-9]+)"), regex.compile(r"(T[0-9]+)"), regex.compile(r"(?<=L)([0-1]+)")))
