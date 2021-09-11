from database_interface import DatabaseInterfaceSequences

database_path = "./results/databases/T1_D80_database_v2.db"
# database_path = snakemake.input[0]

TABLE_ID_SEQ = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    # results = db.query(
    #     f"SELECT id, read_count, sequence_id, sequence, reference_name, cleaved_prefix, ligand_present FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)
    id_sequence_rows = db.get(table=TABLE_ID_SEQ, columns=["id"])

    num_ids = len(id_sequence_rows)
    for rowid in range(num_ids):
        db.update_column_value(TABLE_ID_SEQ, rowid, "possible_sensor", 0)

    # for result in results:
    #     print(result)
