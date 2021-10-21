from database_interface import DatabaseInterfaceSequences
from tqdm import tqdm

database_path = snakemake.input[0]
# database_path = "results/databases/T1_D80_database.db"

# The table where we store all the initial data
TABLE_RAW_SEQ = "raw_sequences"

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"


# Fetch data for raw_sequences
with DatabaseInterfaceSequences(path=database_path) as db:

    # Create the table: 'clean_sequences'

    if not db.table_exists(TABLE_CLEAN_SEQ):
        db.query(f"""CREATE TABLE IF NOT EXISTS {TABLE_CLEAN_SEQ}(
                    id INTEGER PRIMARY KEY,
                    read_count INTEGER,
                    sequence TEXT,
                    sequence_id INTEGER,
                    cleaved_prefix INTEGER,
                    prefix_name TEXT,
                    reference_name TEXT,
                    selection TEXT,
                    driver_round INTEGER,
                    ligand_present INTEGER,
                    cleavage_fraction REAL,
                    fold_change REAL,
                    possible_sensor INTEGER,
                    k_factor REAL,
                    cleavage_fraction_estimated_mean REAL,
                    cleavage_fraction_standard_error REAL,
                    fold_change_estimated_mean REAL,
                    fold_change_standard_error REAL,
                    p_value REAL
                    )""")
    else:
        print(
            f"Table: '{TABLE_CLEAN_SEQ}' already exists, check if files are already processed.")
        while True:
            user_input = input(
                "Do you want to continue filling the database? (Y/n)\t")
            if user_input.lower().startswith('y'):
                print("Database filling up...")
                break
            elif user_input.lower().startswith('n'):
                exit()

    # Retrieve all unique sequences
    id_sequence_rows = db.get(table=TABLE_ID_SEQ, columns=["id", "sequence"])

    for id_sequence_row in tqdm(id_sequence_rows):
        seq_id = id_sequence_row[0]

        for cleaved_bool_val in [True, False]:
            for ligand_bool_val in [True, False]:

                # # Check if read counts were already combined for a sequence id
                # if len(db.retrieve_info_sequence_id(
                #         table=TABLE_CLEAN_SEQ, sequence_id=id_value, cleaved_prefix=cleaved_bool_val, ligand_present=ligand_bool_val)) == 0:

                sequence_rows = db.retrieve_info_sequence_id(
                    table=TABLE_RAW_SEQ, sequence_id=seq_id, cleaved_prefix=cleaved_bool_val, ligand_present=ligand_bool_val)
                if len(sequence_rows) != 0:

                    total_read_count = 0
                    for row in sequence_rows:
                        total_read_count += row[1]

                    sequence_info = {"read_count": total_read_count, "sequence": sequence_rows[0][3], "sequence_id": seq_id, "cleaved_prefix": sequence_rows[0][6], "selection": sequence_rows[0][9], "ligand_present": sequence_rows[0][11],
                                     "prefix_name": sequence_rows[0][7], "reference_name": sequence_rows[0][8], "driver_round": sequence_rows[0][10]}

                    db.query(f"""INSERT INTO {TABLE_CLEAN_SEQ}(read_count, sequence, sequence_id,
                                        cleaved_prefix, prefix_name, reference_name,
                                        selection, driver_round, ligand_present) VALUES (
                                        :read_count, :sequence, :sequence_id,
                                        :cleaved_prefix, :prefix_name, :reference_name,
                                        :selection, :driver_round, :ligand_present)""", parameters=sequence_info)
