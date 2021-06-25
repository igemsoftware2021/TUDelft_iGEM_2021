from pathlib import Path
from tqdm import tqdm
import seq_helpers
import yaml_read_helpers

seq_input_file = snakemake.input[0]
seq_neg_ligand_file, seq_pos_ligand_file = snakemake.output[0], snakemake.output[1]

# Create the path to the config file
config_file_path = Path(__file__).resolve(
).parents[2] / "config" / "config.yaml"

ligand_pos_prefix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
    config_file_path, pattern="prefix", ligand_present=True)
ligand_neg_prefix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
    config_file_path, pattern="prefix", ligand_present=False)

with open(seq_input_file, 'r') as rf:
    lines = rf.readlines()
    with open(seq_pos_ligand_file, "w") as pf, open(seq_neg_ligand_file, "w") as nf:
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
