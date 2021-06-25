from database_interface import DatabaseInterfaceCleanSequences
from tqdm import tqdm

# database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
database_path = snakemake.input[0]

# Fetch data for raw_sequences
with DatabaseInterfaceCleanSequences(path=database_path) as db:
    # Fetchall rows only columns (cleaned_sequence, prefix_name, ligand_present)
    pre_data_seq = db.get(table="raw_sequences", columns=[
                          "cleaned_sequence", "cleaved_prefix", "ligand_present"])

    # Create a set of the list of tuples of the retrieved sequence data
    unique_rows = list(set(pre_data_seq))

TABLE_NAME = "clean_sequences"

# Create table
with DatabaseInterfaceCleanSequences(path=database_path) as db:

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

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    for i in tqdm(range(len(unique_rows))):
        # get for every unique sequence all the read_counts
        sequence_all = db.get_info_sequence(
            "raw_sequences", unique_rows[i][0], unique_rows[i][1], unique_rows[i][2])

        read_counts_new = sum([p[1] for p in sequence_all])

        clean_sequence_info = {"read_count": read_counts_new, "cleaned_sequence": sequence_all[0][3], "cleaved_prefix": sequence_all[0][5], "selection": sequence_all[0][8], "ligand_present": sequence_all[0][10],
                               "prefix_name": sequence_all[0][6], "reference_name": sequence_all[0][7], "driver_round": sequence_all[0][9], "cleavage_fraction": "NULL", "fold_change": "NULL", "possible_sensor": 0,
                               "k_factor": "NULL", "cleavage_fraction_estimated_mean": "NULL", "cleavage_fraction_standard_deviation": "NULL", "fold_change_estimated_mean": "NULL", "fold_change_standard_deviation": "NULL",
                               "fold_change_standard_error": "NULL"}

        # insert the information into the new table
        db.insert_sequence_info(TABLE_NAME, clean_sequence_info)
