from pathlib import Path
import yaml
import regex
from tqdm import tqdm
from database_interface import DatabaseInterfaceSequences
import seq_helpers
import yaml_read_helpers

# Find the path to the config file
config_file_path = Path(__file__).resolve(
).parents[2] / "config" / "config_folate.yaml"

# Find the path to the references file
ngs_references_path = Path(
    __file__).resolve().parents[2] / "data" / "ngs_references.csv"

# Store the input filenames in a variable. The inputfile names are given by snakemake.
# inputfiles = [
# Path(__file__).resolve().parents[2] / "results" /
# "read_counts" / "T1_D80_L0_read_count.txt",
# Path(__file__).resolve().parents[2] / "results" / "read_counts" / "T1_D80_L1_read_count.txt"]
inputfiles = [snakemake.input[0], snakemake.input[1]]

# Store the output filename in a variable. The outputfile name is given by snakemake.
# database_path = ":memory:"
# database_path = Path(__file__).resolve(
# ).parents[2] / "results" / "databases" / "T1_D80_database.db"
database_path = snakemake.output[0]

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


print(clean_ngs_reference_patterns)

# Retrieve the compiled info patterns from the config file
driver_round_pattern, selection_pattern, ligand_present_pattern = yaml_read_helpers.retrieve_compiled_info_patterns(
    config_file_path)

minimum_number_reads = yaml_read_helpers.retrieve_minimum_reads(
    config_file_path)

# The table where we store all the initial data
TABLE_RAW_SEQ = "raw_sequences"

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

print(database_path)
# Create the database
with DatabaseInterfaceSequences(path=database_path) as db:

    if not db.table_exists(TABLE_RAW_SEQ):
        db.query(f"""CREATE TABLE IF NOT EXISTS {TABLE_RAW_SEQ}(
                    id INTEGER PRIMARY KEY,
                    read_count INTEGER,
                    original_sequence TEXT,
                    sequence TEXT,
                    sequence_id INTEGER,
                    barcode TEXT,
                    cleaved_prefix INTEGER,
                    prefix_name TEXT,
                    reference_name TEXT,
                    selection TEXT,
                    driver_round INTEGER,
                    ligand_present INTEGER
                    )""")
    else:
        print(
            f"Table: '{TABLE_RAW_SEQ}' already exists, check if files are already processed.")
        while True:
            user_input = input(
                "Do you want to continue filling the database? (Y/n)\t")
            if user_input.lower().startswith('y'):
                print("Database filling up...")
                break
            elif user_input.lower().startswith('n'):
                exit()

    if not db.table_exists(TABLE_ID_SEQ):
        db.query(f"""CREATE TABLE IF NOT EXISTS {TABLE_ID_SEQ}(
                    id INTEGER PRIMARY KEY,
                    sequence TEXT,
                    reference_name TEXT,
                    negative_cs REAL,
                    positive_cs REAL,
                    fold_change REAL
                    )""")
    else:
        print(
            f"Table: '{TABLE_ID_SEQ}' already exists, check if files are already processed.")
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
        inputfile = str(inputfile)
        print(f"Now inserting the sequences from {inputfile}")

        # Retrieve the round of DRIVER, which selection it is and whether the ligand
        # is present from the file
        print(driver_round_pattern, inputfile)
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

                read_count, sequence = line.strip().split()
                if int(read_count) >= minimum_number_reads:

                    # Create a dictionary and store all general information for an unique sequence
                    sequence_info = {"driver_round": driver_round,
                                     "selection": selection, "ligand_present": ligand_present}
                    # sequence_info["reference_name"] = seq_helpers.reference_seq(
                    #     sequence, clean_ngs_reference_patterns)

                    # if sequence_info["reference_name"] is not None:
                    #     print('yeahhhh', sequence_info["reference_name"])

                    sequence_info["read_count"] = int(read_count)
                    sequence_info["original_sequence"] = sequence

                    # Determine the prefix sequence and the corresponding name
                    prefix_seq, prefix_name = seq_helpers.determine_pattern(
                        sequence, patterns_info=prefix_patterns)

                    # Determine the suffix sequence and the corresponding name
                    suffix_seq, suffix_name = seq_helpers.determine_pattern(
                        sequence, patterns_info=suffix_patterns)

                    clvd_prefix = seq_helpers.determine_clvd_prefix(
                        prefix_name, clvd_prefix_name=clvd_prefix_name, unclvd_prefix_name=unclvd_prefix_name)

                    sequence_info["cleaved_prefix"] = clvd_prefix
                    sequence_info["prefix_name"] = prefix_name

                    sequence_info["barcode"] = seq_helpers.retrieve_barcode(
                        sequence, prefix_seq)

                    sequence_info["sequence"] = seq_helpers.clean_sequence(
                        sequence, prefix_seq, suffix_seq)

                    if sequence_info["sequence"] is not None:
                        # Determine whether the cleaned sequence is a cleaned reference sequence
                        sequence_info["reference_name"] = seq_helpers.reference_seq(
                            sequence_info["sequence"], clean_ngs_reference_patterns)

                        if sequence_info["reference_name"] is not None:
                            print('yeahhhh', sequence_info["reference_name"])

                        # This part is to determine the unique sequence id for a certain sequence
                        id_sequence_info = db.retrieve_info_sequence(
                            table=TABLE_ID_SEQ, sequence=sequence_info["sequence"])

                        # print(id_sequence_info)

                        if len(id_sequence_info) == 0:
                            db.query(f"""INSERT INTO {TABLE_ID_SEQ} (sequence) VALUES (:sequence)""", parameters={
                                "sequence": sequence_info["sequence"]})
                            sequence_id = db.cursor.lastrowid  # rowid in sequence table

                            if sequence_info["reference_name"] is not None:
                                db.update_column_value(
                                    TABLE_ID_SEQ, rowid=sequence_id, column_name="reference_name", value=sequence_info["reference_name"])
                        else:
                            sequence_id = id_sequence_info[0][0]

                        sequence_info["sequence_id"] = sequence_id

                        db.query(f"""INSERT INTO {TABLE_RAW_SEQ} (read_count, original_sequence, sequence, sequence_id,
                                barcode, cleaved_prefix, prefix_name, reference_name,
                                selection, driver_round, ligand_present) VALUES (
                                :read_count, :original_sequence, :sequence, :sequence_id,
                                :barcode, :cleaved_prefix, :prefix_name, :reference_name,
                                :selection, :driver_round, :ligand_present)""", parameters=sequence_info)


# with DatabaseInterfaceRawSequences(path=database_path) as db:
#     results = db.get(
#         TABLE_NAME, ["original_sequence", "cleaned_sequence"], limit=10)

#     # testing = db.get_ref_sequences()
#     # print(testing)
#     # print(db.get_sequences(ligand_present=0))
