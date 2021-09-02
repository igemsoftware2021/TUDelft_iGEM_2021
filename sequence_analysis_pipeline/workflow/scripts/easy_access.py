from database_interface import DatabaseInterfaceCleanSequences

database_path = "./results/databases/T1_D80_database.db"
# database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    sequences_info = db.query(
        f"SELECT cleaned_sequence FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)
    uniq_sequences_info = list(set(sequences_info))

print(len(uniq_sequences_info))
