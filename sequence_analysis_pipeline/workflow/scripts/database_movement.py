from database_interface import DatabaseInterfaceCleanSequences


database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
# database_path = snakemake.output[0]

# Fetch data for raw_sequences
with DatabaseInterface(path=database_path) as db:
    # Fetchall rows only columns (cleaned_sequence, prefix_name, ligand_present)
    pre_data_seq = db.get(table="raw_sequences", columns=[
                          "cleaned_sequence", "cleaved_prefix", "ligand_present"])

    # Create a set of the list of tuples of the retrieved sequence data
    unique_rows = list(set(pre_data_seq))


# Create table
with DatabaseInterfaceCleanSequences(path=database_path) as db:
    TABLE_NAME = "clean_sequences"

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
    for i in range(len(unique_rows)):
        # get for every unique sequence all the read_counts
        sequence_all = db.get_all_unique_sequence(table="raw_sequences", cleaned_sequences=unique_rows[i][0], cleaved_prefix=unique_rows[i][1], ligand_present=unique_rows[i][2]])
        read_counts_new= sum([p[1] for p in sequence_all])

# Go over all the tuples in the set, query the data, add the read counts and put all info
# in new table.
# Then do the cleavage fraction and stuff
