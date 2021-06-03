import yaml
import regex
from tqdm import tqdm
from database_interface import DatabaseInterfaceSequences
import seq_helper


# Read out all the information needed from the 'config.yaml' file
with open("sequence_analysis_pipeline/config.yaml", "r") as rf:
    try:
        yaml_info = yaml.safe_load(rf)
        clvd_prefix_seq = yaml_info["prefix"]["cleaved"]["sequence"]
        clvd_prefix_name = yaml_info["prefix"]["cleaved"]["name"]
        clvd_suffix_seq = yaml_info["suffix"]["cleaved"]["sequence"]

        unclvd_prefix_seq = yaml_info["prefix"]["uncleaved"]["sequence"]
        unclvd_prefix_name = yaml_info["prefix"]["uncleaved"]["name"]
        unclvd_suffix_seq = yaml_info["suffix"]["uncleaved"]["sequence"]

        # Store the clvd and unclvd prefix info in a dictionary for easy reference
        clvd_prefix_info = {"seq": clvd_prefix_seq, "name": clvd_prefix_name}
        unclvd_prefix_info = {
            "seq": unclvd_prefix_seq, "name": unclvd_prefix_name}

        driver_round_pattern = yaml_info["info_patterns"]["driver_round"]
        selection_pattern = yaml_info["info_patterns"]["selection"]
        ligand_present_pattern = yaml_info["info_patterns"]["ligand_present"]

    except yaml.YAMLError as exc:
        print(exc)

# Filename is given by snakemake
inputfiles = [
    'sequence_analysis_pipeline/data/NGS/processed/S1_D80_L0_read_count.txt']
#inputfiles = [snakemake.input[0], snakemake.input[1]]

# database_path = ":memory:"
database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
# database_path = snakemake.output[0]

ngs_references = seq_helper.read_ngs_references(
    "sequence_analysis_pipeline/ngs_references.csv")

# Complement, reverse and cleanup the ngs reference sequences. Cleaning up
# means removing the prefix and suffix from the sequence.
ready_ngs_references = seq_helper.clean_ngs_reference_sequences(
    ngs_references)

# Create a dictionary with the sequences as compiled regex objects. This is
# for optimization.
clean_ngs_reference_patterns = seq_helper.create_ngs_references_patterns(
    ready_ngs_references)

# Create the database
with DatabaseInterfaceSequences(path=database_path) as db:

    if not db.table_exists("raw_sequences"):
        db.create_table("raw_sequences")
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


with DatabaseInterfaceSequences(path=database_path) as db:
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
                                 "cleavage_fraction": "NULL", "fold_change": "NULL", "possible_sensor": 0, "mutated_prefix": 0}

                read_count, sequence = line.strip().split()
                sequence_info["read_count"] = read_count
                sequence_info["original_sequence"] = sequence

                prefix_seq, prefix_name = seq_helper.determine_prefix(sequence)

                clvd_prefix = seq_helper.determine_clvd_prefix(
                    prefix_name, clvd_prefix_name=clvd_prefix_name, unclvd_prefix_name=unclvd_prefix_name)

                sequence_info["cleaved_prefix"] = clvd_prefix
                sequence_info["prefix_name"] = prefix_name

                sequence_info["barcode"] = seq_helper.retrieve_barcode(
                    sequence, prefix_seq)

                sequence_info["cleaned_sequence"] = seq_helper.clean_sequence(
                    sequence, prefix_seq, clvd_suffix_seq)

                # Determine whether the cleaned sequence is a cleaned reference sequence
                sequence_info["reference_name"] = seq_helper.reference_seq(
                    sequence_info["cleaned_sequence"], clean_ngs_reference_patterns)

                db.insert_sequence_info("raw_sequences", sequence_info)

with DatabaseInterfaceSequences(path=database_path) as db:
    results = db.get(
        "sequences", ["original_sequence", "cleaned_sequence"], limit=10)

    print(results)
    # testing = db.get_ref_sequences()
    # print(testing)
    # print(db.get_sequences(ligand_present=0))


