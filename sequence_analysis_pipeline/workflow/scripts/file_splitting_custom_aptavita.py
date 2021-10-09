from pathlib import Path
from tqdm import tqdm
import seq_helpers
import yaml_read_helpers

# seq_input_file = "results/read_counts/N35-I3_S3.txt"
seq_input_file = snakemake.input[0]
# seq_neg_ligand_file, seq_pos_ligand_file = "results/read_counts/N35-I3_S3_L0_read_count.txt", "results/read_counts/N35-I3_S3_L0_read_count.txt"
seq_neg_ligand_file, seq_pos_ligand_file, seq_other_file = snakemake.output[
    0], snakemake.output[1], snakemake.output[2]

# Create the path to the config file
config_file_path = Path(__file__).resolve(
).parents[2] / "config" / "config.yaml"

ligand_pos_prefix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
    config_file_path, pattern="prefix", ligand_present=True)
ligand_pos_suffix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
    config_file_path, pattern="suffix", ligand_present=True)

ligand_neg_prefix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
    config_file_path, pattern="prefix", ligand_present=False)
ligand_neg_suffix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
    config_file_path, pattern="suffix", ligand_present=False)

with open(seq_input_file, 'r') as rf:
    lines = rf.readlines()
    with open(seq_pos_ligand_file, "w") as pf, open(seq_neg_ligand_file, "w") as nf, open(seq_other_file, "w") as of:
        for line in tqdm(lines):
            read_count, sequence = line.strip().split()
            prefix_seq, pattern_name_prefix = seq_helpers.determine_pattern(
                sequence, ligand_pos_prefix_patterns)
            suffix_seq, pattern_name_suffix = seq_helpers.determine_pattern(
                sequence, ligand_pos_suffix_patterns)

            if pattern_name_prefix is not None and pattern_name_suffix is not None:
                pf.write(read_count + " " + sequence + "\n")
            else:
                prefix_seq, pattern_name_prefix = seq_helpers.determine_pattern(
                    sequence, ligand_neg_prefix_patterns)
                suffix_seq, pattern_name_suffix = seq_helpers.determine_pattern(
                    sequence, ligand_neg_suffix_patterns)
                if pattern_name_prefix is not None and pattern_name_suffix is not None:
                    nf.write(read_count + " " + sequence + "\n")

                # If it does not belong to either positive or negative ligand
                # store it in the "other"-file.
                else:
                    of.write(read_count + " " + sequence + "\n")
