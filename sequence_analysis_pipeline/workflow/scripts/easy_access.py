from database_interface import DatabaseInterfaceSequences

database_path = "./results/databases/T1_D80_database.db"
# database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    results = db.query(
        f"SELECT id, read_count, sequence_id, sequence, cleaved_prefix, ligand_present FROM {TABLE_NAME} WHERE reference_name IS NOT NULL", fetchall=True)

    for result in results:
        print(result)
