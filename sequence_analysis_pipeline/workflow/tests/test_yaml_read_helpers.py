import unittest
import regex
from workflow.scripts.yaml_read_helpers import retrieve_dataset_name
from workflow.scripts.yaml_read_helpers import retrieve_minimum_reads
from workflow.scripts.yaml_read_helpers import retrieve_compiled_patterns
from workflow.scripts.yaml_read_helpers import retrieve_compiled_reference_patterns
from workflow.scripts.yaml_read_helpers import retrieve_prefix_name
from workflow.scripts.yaml_read_helpers import retrieve_compiled_info_patterns


class TestSuiteCalcHelpers(unittest.TestCase):
    """Test suite to test the calc helpers functions"""

    def test_retrieve_dataset_name(self):

        yaml_file_path = "./config/config_example.yaml"

        self.assertEqual(retrieve_dataset_name(yaml_file_path), "S1_D63")

    def test_retrieve_minimum_reads(self):
        yaml_file_path = "./config/config_example.yaml"

        self.assertEqual(retrieve_minimum_reads(yaml_file_path), 25)

    def test_retrieve_compiled_patterns(self):
        """Function tests the retrieve_compiled_patterns() function."""

        yaml_file_path = "./config/config_example.yaml"

        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="prefix", ligand_present=True),
            {regex.compile("(?e)(TATCATGCAGAATTTAATACGACTCACTATAGGGACAAAACAAAACG){e<=1}"): "Z2",
             regex.compile("(?e)(GATCTCATTCAATTTAATACGACTCACTATAGGGAAACAAACAAAG){e<=1}"): "W2"})
        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="suffix", ligand_present=True),
            {regex.compile("(?e)(AAAAAGAAAAATAAAAAGTACGCAGATC){e<=1}"): "X1"})
        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="prefix", ligand_present=False),
            {regex.compile("(?e)(CGCTTATCCTAATTTAATACGACTCACTATAGGGACAAAACAAAACG){e<=1}"): "Z1",
             regex.compile("(?e)(AGTCATTGAGAATTTAATACGACTCACTATAGGGAAACAAACAAAG){e<=1}"): "W1"})
        self.assertEqual(retrieve_compiled_patterns(
            yaml_file_path, pattern="suffix", ligand_present=False),
            {regex.compile("(?e)(AAAAAGAAAAATAAAAAGTACGCAGATC){e<=1}"): "X1"})

    def test_retrieve_compiled_reference_patterns(self):
        """Function tests the retrieve_compiled_patterns() function."""

        yaml_file_path = "./config/config_example.yaml"

        self.assertEqual(retrieve_compiled_reference_patterns(
            yaml_file_path, pattern="prefix"),
            {regex.compile("(?e)(ACAAAACAAAAC){e<=1}"): "Z",
             regex.compile("(?e)(GGGAAACAAACAAA){e<=1}"): "W"})
        self.assertEqual(retrieve_compiled_reference_patterns(
            yaml_file_path, pattern="suffix"),
            {regex.compile("(?e)(AAAAAGAAA){e<=1}"): "X"})

    def test_retrieve_prefix_name(self):
        """Function tests the retrieve_prefix_name() function."""

        yaml_file_path = "./config/config_example.yaml"

        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=True, ligand_present=True), "Z2")
        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=False, ligand_present=True), "W2")
        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=True, ligand_present=False), "Z1")
        self.assertEqual(retrieve_prefix_name(
            yaml_file_path, cleaved=False, ligand_present=False), "W1")

    def test_retrieve_compiled_info_patterns(self):
        """Function tests the retrieve_compiled_info_patterns() function."""

        yaml_file_path = "./config/config_example.yaml"

        self.assertEqual(retrieve_compiled_info_patterns(yaml_file_path), (regex.compile(
            r"(?<=S)([0-9]+)"), regex.compile(r"(D[0-9]+)"), regex.compile(r"(?<=L)([0-1]+)")))
