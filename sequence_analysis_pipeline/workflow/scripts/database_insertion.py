from pathlib import Path
import yaml
import regex
from tqdm import tqdm
from database_interface import DatabaseInterfaceRawSequences
import seq_helpers
import yaml_read_helpers

# Find the path to the config file
config_file_path = Path(__file__).resolve(
).parents[2] / "config" / "config.yaml"

# Find the path to the references file
ngs_references_path = config_file_path = Path(
    __file__).resolve().parents[2] / "data" / "ngs_references.csv"

# Store the input filenames in a variable. The inputfile names are given by snakemake.
inputfiles = [
    'sequence_analysis_pipeline/data/NGS/T1_D80_L0_read_count.txt']
# inputfiles = [snakemake.input[0], snakemake.input[1]]

# Store the output filename in a variable. The outputfile name is given by snakemake.
# database_path = ":memory:"
database_path = "sequence_analysis_pipeline/data/NGS/T1_D80_database.db"
# database_path = snakemake.output[0]

# Prep the reference sequences
reference_prefix_patterns = yaml_read_helpers.retrieve_compiled_reference_patterns(
    config_file_path, pattern="prefix")
reference_suffix_patterns = yaml_read_helpers.retrieve_compiled_reference_patterns(
    config_file_path, pattern="suffix")

ngs_references = seq_helpers.read_ngs_references(ngs_references_path)

# Complement, reverse and cleanup the ngs reference sequences. Cleaning up
# means removing the prefix and suffix from the sequence.
ready_ngs_references = seq_helpers.clean_ngs_reference_sequences(
    ngs_references, reference_prefix_patterns, reference_suffix_patterns)

# Create a dictionary with the sequences as compiled regex objects. This is
# for optimization.
clean_ngs_reference_patterns = seq_helpers.create_ngs_references_patterns(
    ready_ngs_references)

# Retrieve the compiled info patterns from the config file
driver_round_pattern, selection_pattern, ligand_present_pattern = yaml_read_helpers.retrieve_compiled_info_patterns(
    config_file_path)

# The name of the database table is a constant, because it is the same for all the operations.
TABLE_NAME = "raw_sequences"

# Create the database
with DatabaseInterfaceRawSequences(path=database_path) as db:

    if not db.table_exists(TABLE_NAME):
        db.create_table(TABLE_NAME)
    else:
        print("Table already exists, check if files are already processed.")
        while True:
            user_input = input(
                "Do you want to continue filling the database? (Y/n)\t")
            if user_input.lower().startswith('y'):
                print("Database filling up...")
                break
            elif user_input.lower().startswith('n'):
                exit()

with DatabaseInterfaceRawSequences(path=database_path) as db:
    print("Insertion of the sequences in the database...")
    for inputfile in inputfiles:
        print(f"Now inserting the sequences from {inputfile}")

        # Retrieve the round of DRIVER, which selection it is and whether the ligand
        # is present from the file
        driver_round = int(driver_round_pattern.search(inputfile).group())
        selection = selection_pattern.search(inputfile).group()
        ligand_present = int(ligand_present_pattern.search(inputfile).group())

        prefix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
            config_file_path, pattern="prefix", ligand_present=bool(ligand_present))
        suffix_patterns = yaml_read_helpers.retrieve_compiled_patterns(
            config_file_path, pattern="suffix", ligand_present=bool(ligand_present))

        clvd_prefix_name = yaml_read_helpers.retrieve_prefix_name(
            config_file_path, cleaved=True, ligand_present=bool(ligand_present))
        unclvd_prefix_name = yaml_read_helpers.retrieve_prefix_name(
            config_file_path, cleaved=False, ligand_present=bool(ligand_present))

        with open(inputfile) as rf:
            lines = rf.readlines()
            for line in tqdm(lines):
                # Create a dictionary and store all general information for an unique sequence
                sequence_info = {"driver_round": driver_round, "selection": selection, "ligand_present": ligand_present,
                                 "cleavage_fraction": "NULL", "fold_change": "NULL", "possible_sensor": 0}

                read_count, sequence = line.strip().split()
                sequence_info["read_count"] = int(read_count)
                sequence_info["original_sequence"] = sequence

                # Determine the prefix sequence and the corresponding name
                prefix_seq, prefix_name, mutated_prefix = seq_helpers.determine_pattern(
                    sequence, patterns_info=prefix_patterns)
                sequence_info["mutated_prefix"] = mutated_prefix

                # Determine the suffix sequence and the corresponding name
                suffix_seq, suffix_name, mutated_suffix = seq_helpers.determine_pattern(
                    sequence, patterns_info=suffix_patterns)
                sequence_info["mutated_suffix"] = mutated_suffix

                clvd_prefix = seq_helpers.determine_clvd_prefix(
                    prefix_name, clvd_prefix_name=clvd_prefix_name, unclvd_prefix_name=unclvd_prefix_name)

                sequence_info["cleaved_prefix"] = clvd_prefix
                sequence_info["prefix_name"] = prefix_name

                sequence_info["barcode"] = seq_helpers.retrieve_barcode(
                    sequence, prefix_seq)

                sequence_info["cleaned_sequence"] = seq_helpers.clean_sequence(
                    sequence, prefix_seq, suffix_seq)

                # Determine whether the cleaned sequence is a cleaned reference sequence
                sequence_info["reference_name"] = seq_helpers.reference_seq(
                    sequence_info["cleaned_sequence"], clean_ngs_reference_patterns)

                db.insert_sequence_info(TABLE_NAME, sequence_info)


# with DatabaseInterfaceRawSequences(path=database_path) as db:
#     results = db.get(
#         TABLE_NAME, ["original_sequence", "cleaned_sequence"], limit=10)

#     # testing = db.get_ref_sequences()
#     # print(testing)
#     # print(db.get_sequences(ligand_present=0))
