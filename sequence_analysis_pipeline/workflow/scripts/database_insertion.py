from pathlib import Path
import yaml
import regex
from tqdm import tqdm
from database_interface import DatabaseInterfaceRawSequences
import seq_helper

# Create the path to the config file
config_file_path = Path(__file__).resolve(
).parents[2] / "config" / "config.yaml"

# Read out all the information needed from the 'config.yaml' file
with open(config_file_path, "r") as rf:
    try:
        yaml_info = yaml.safe_load(rf)
        prefixes = yaml_info["prefixes"]

        # Create regex patterns from the prefix sequences
        prefix_patterns = {}
        for prefix in prefixes:
            sequence = yaml_info["prefixes"][prefix]["sequence"]
            prefix_name = yaml_info["prefixes"][prefix]["name"]
            max_error = yaml_info["prefixes"][prefix]["max_error"]
            prefix_patterns[regex.compile(
                fr"(?e)({sequence}){{e<={max_error}}}")] = prefix_name

        # Create a dictionary with regex patterns for the suffix sequences
        suffixes = yaml_info["suffixes"]
        suffix_patterns = {}
        for suffix in suffixes:
            sequence = yaml_info["suffixes"][suffix]["sequence"]
            suffix_name = yaml_info["suffixes"][suffix]["name"]
            max_error = yaml_info["suffixes"][suffix]["max_error"]
            suffix_patterns[regex.compile(
                fr"(?e)({sequence}){{e<={max_error}}}")] = suffix_name

        # Store the name of the cleaved prefix
        clvd_prefix_name = yaml_info["prefixes"]["cleaved"]["name"]

        # Store the name of the uncleaved prefix
        unclvd_prefix_name = yaml_info["prefixes"]["uncleaved"]["name"]

        driver_round_pattern = yaml_info["info_patterns"]["driver_round"]
        selection_pattern = yaml_info["info_patterns"]["selection"]
        ligand_present_pattern = yaml_info["info_patterns"]["ligand_present"]

    except yaml.YAMLError as exc:
        print(exc)

# Filename is given by snakemake
inputfiles = [
    'sequence_analysis_pipeline/data/NGS/T1_D80_L0_read_count.txt']
# inputfiles = [snakemake.input[0], snakemake.input[1]]

# database_path = ":memory:"
database_path = "sequence_analysis_pipeline/data/NGS/T1_D80_database.db"
# database_path = snakemake.output[0]

ngs_references_path = config_file_path = Path(
    __file__).resolve().parents[2] / "data" / "ngs_references.csv"

ngs_references = seq_helper.read_ngs_references(ngs_references_path)

# Complement, reverse and cleanup the ngs reference sequences. Cleaning up
# means removing the prefix and suffix from the sequence.
ready_ngs_references = seq_helper.clean_ngs_reference_sequences(
    ngs_references, prefix_patterns, suffix_patterns)

# Create a dictionary with the sequences as compiled regex objects. This is
# for optimization.
clean_ngs_reference_patterns = seq_helper.create_ngs_references_patterns(
    ready_ngs_references)

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
        driver_round = int(regex.search(
            driver_round_pattern, inputfile).group())
        selection = regex.search(selection_pattern, inputfile).group()
        ligand_present = int(
            regex.search(ligand_present_pattern, inputfile).group())

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
                prefix_seq, prefix_name, mutated_prefix = seq_helper.determine_pattern(
                    sequence, patterns_info=prefix_patterns)
                sequence_info["mutated_prefix"] = mutated_prefix

                # Determine the suffix sequence and the corresponding name
                suffix_seq, suffix_name, mutated_suffix = seq_helper.determine_pattern(
                    sequence, patterns_info=suffix_patterns)
                sequence_info["mutated_suffix"] = mutated_suffix

                clvd_prefix = seq_helper.determine_clvd_prefix(
                    prefix_name, clvd_prefix_name=clvd_prefix_name, unclvd_prefix_name=unclvd_prefix_name)

                sequence_info["cleaved_prefix"] = clvd_prefix
                sequence_info["prefix_name"] = prefix_name

                sequence_info["barcode"] = seq_helper.retrieve_barcode(
                    sequence, prefix_seq)

                sequence_info["cleaned_sequence"] = seq_helper.clean_sequence(
                    sequence, prefix_seq, suffix_seq)

                # Determine whether the cleaned sequence is a cleaned reference sequence
                sequence_info["reference_name"] = seq_helper.reference_seq(
                    sequence_info["cleaned_sequence"], clean_ngs_reference_patterns)

                db.insert_sequence_info(TABLE_NAME, sequence_info)


with DatabaseInterfaceRawSequences(path=database_path) as db:
    results = db.get(
        TABLE_NAME, ["original_sequence", "cleaned_sequence"], limit=10)

    # testing = db.get_ref_sequences()
    # print(testing)
    # print(db.get_sequences(ligand_present=0))
