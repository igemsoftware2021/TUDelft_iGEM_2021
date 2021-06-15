import numpy as np
from database_interface import DatabaseInterfaceCleanSequences
import bootstrapping

database_path = "sequence_analysis_pipeline/data/NGS/T1_D80_database.db"
# database_path = snakemake.output[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    ref_unclvd_neg = db.get_ref_sequences(
        table=TABLE_NAME, cleaved_prefix=0, ligand_present=0)
    r_ref_unclvd_neg = sum(seq_info[1] for seq_info in ref_unclvd_neg)

    # First retrieve all sequences
    sequences_all = db.get(table=TABLE_NAME, columns=["cleaned_sequence"])
    # Remove duplicate sequences
    sequences_unique = list(set(sequences_all))

    for sequence in sequences_unique:
        # r_ = read_count
        r_seq_unclvd_neg = db.get_info_sequence(
            table=TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=0, ligand_present=0)[1]
        r_seq_clvd_neg = db.get_info_sequence(
            table=TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=1, ligand_present=0)[1]

        data = np.zeros(r_seq_unclvd_neg+r_seq_clvd_neg, dtype=np.int8)
        data[-r_seq_clvd_neg:] = 1
