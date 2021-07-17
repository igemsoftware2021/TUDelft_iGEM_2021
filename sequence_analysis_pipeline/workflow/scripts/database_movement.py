from database_interface import DatabaseInterfaceCleanSequences
from tqdm import tqdm

# database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
database_path = snakemake.input[0]
# database_path = "results/databases/T1_D80_database.db"

# Store in this variable all tuples of ("cleaned_sequence", "cleaved_prefix", "ligand_present") that have already been checked
sequence_info_checked = set()

# Fetch data for raw_sequences
with DatabaseInterfaceCleanSequences(path=database_path) as db:
    # Fetchall all rows only columns (id, cleaned_sequence, prefix_name, ligand_present)
    rows_info_seq = db.get(table="raw_sequences", columns=[
        "id", "cleaned_sequence", "cleaved_prefix", "ligand_present"])

    # Fetchall all rows only columns (id, read_count)
    rows_id_and_read_count = db.get(table="raw_sequences", columns=[
        "id", "read_count"])

# Create dictionary where (sequence, cleaved_prefix, ligand_present) is the key
# and the value is a list of all the rowids corresponding to this
info_seq_rowid_dict = {}
for row_info in rows_info_seq:
    key_tuple = (row_info[1], row_info[2], row_info[3])
    rowids = info_seq_rowid_dict.get(key_tuple, [])
    rowids.append(row_info[0])
    info_seq_rowid_dict[key_tuple] = rowids

# Dictionary where key is the rowid, and the value is the read count
id_read_count_dict = {row_info[0]: row_info[1]
                      for row_info in rows_id_and_read_count}


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

# This has a speed of approximately 9000 sequences per second
with DatabaseInterfaceCleanSequences(path=database_path) as db:
    for key in tqdm(info_seq_rowid_dict):
        total_read_count = 0
        rowids = info_seq_rowid_dict[key]
        for rowid in rowids:
            total_read_count += id_read_count_dict[rowid]

        # Retrieve the sequence info of a certain sequence
        sequence_info = db.query(
            f"SELECT * FROM raw_sequences WHERE id={rowids[0]}", fetchall=True)

        clean_sequence_info = {"read_count": total_read_count, "cleaned_sequence": sequence_info[0][3], "cleaved_prefix": sequence_info[0][5], "selection": sequence_info[0][8], "ligand_present": sequence_info[0][10],
                               "prefix_name": sequence_info[0][6], "reference_name": sequence_info[0][7], "driver_round": sequence_info[0][9]}

        # insert the information into the new table
        db.insert_movement_sequence_info(TABLE_NAME, clean_sequence_info)


# OLD: this did work but it only accomplished a speed of 6 sequences a second, which is waaaay too slow.
# with DatabaseInterfaceCleanSequences(path=database_path) as db:
#     for i in tqdm(range(len(unique_rows))):
#         # get for every unique sequence all the read_counts
#         sequence_all = db.get_info_sequence(
#             "raw_sequences", unique_rows[i][0], unique_rows[i][1], unique_rows[i][2])

#         read_counts_new = sum([p[1] for p in sequence_all])

#         clean_sequence_info = {"read_count": read_counts_new, "cleaned_sequence": sequence_all[0][3], "cleaved_prefix": sequence_all[0][5], "selection": sequence_all[0][8], "ligand_present": sequence_all[0][10],
#                                "prefix_name": sequence_all[0][6], "reference_name": sequence_all[0][7], "driver_round": sequence_all[0][9], "cleavage_fraction": "NULL", "fold_change": "NULL", "possible_sensor": 0,
#                                "k_factor": "NULL", "cleavage_fraction_estimated_mean": "NULL", "cleavage_fraction_standard_deviation": "NULL", "fold_change_estimated_mean": "NULL", "fold_change_standard_deviation": "NULL",
#                                "fold_change_standard_error": "NULL"}

#         # insert the information into the new table
#         db.insert_sequence_info(TABLE_NAME, clean_sequence_info)
