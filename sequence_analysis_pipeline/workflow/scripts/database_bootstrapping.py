from database_interface import DatabaseInterfaceCleanSequences
import bootstrapping

database_path = "sequence_analysis_pipeline/data/NGS/T1_D80_database.db"
# database_path = snakemake.output[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:

    # First retrieve all sequences
    sequences_all = db.get(table=TABLE_NAME, columns=["clean_sequence"])
    # Remove duplicate sequences
    sequences_unique = list(set(sequences_all))

    for sequence in sequences_unique:
