from database_interface import DatabaseInterfaceCleanSequences

database_path = "results/databases/T1_D80_database.db"

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    cs = db.query(
        f"SELECT cleavage_fraction FROM {TABLE_NAME} WHERE cleavage_fraction IS NOT NULL", fetchall=True)

    fc = db.query(
        f"SELECT fold_change FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)

cs = [info[0] for info in cs]
fc = [info[0] for info in fc]

print(set(fc))
# print(cs)
