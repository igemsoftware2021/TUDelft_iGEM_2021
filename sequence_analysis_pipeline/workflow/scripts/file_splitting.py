import yaml
from pathlib import Path
import regex
from tqdm import tqdm
import seq_helpers
import yaml_reading

# Create the path to the config file
config_file_path = Path(__file__).resolve(
).parents[2] / "config" / "config_splitting.yaml"

ligand_pos_prefix_patterns, ligand_pos_suffix_patterns = yaml_reading.retrieve_compiled_prefix_patterns(
    config_file_path, ligand_present=True)
ligand_neg_prefix_patterns, ligand_neg_suffix_patterns = yaml_reading.retrieve_compiled_prefix_patterns(
    config_file_path, ligand_present=False)


raw_data_file = "data/NGS/N35-I3_S3_read_count.txt"

with open(raw_data_file, 'r') as rf:
    lines = rf.readlines()
    with open("data/NGS/T2_D80_L1.txt", "w") as pf, open("data/NGS/T2_D80_L0.txt", "w") as nf:
        for line in tqdm(lines):
            read_count, sequence = line.strip().split()
            pattern_seq, pattern_name, mutated_pattern = seq_helpers.determine_pattern(
                sequence, ligand_pos_prefix_patterns)
            if pattern_name is not None:
                pf.write(read_count + " " + sequence + "\n")
            else:
                pattern_seq, pattern_name, mutated_pattern = seq_helpers.determine_pattern(
                    sequence, ligand_neg_prefix_patterns)
                if pattern_name is not None:
                    nf.write(read_count + " " + sequence + "\n")