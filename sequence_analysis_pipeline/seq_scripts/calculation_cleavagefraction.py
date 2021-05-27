import sqlite3
from database_interface import DatabaseInterfaceSequences

database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
# database_path = snakemake.input[0]

with DatabaseInterfaceSequences(path=database_path) as db:
    # 1 - Get number of reads of the references with cleaved prefix
    clv_ref_info = db.get_ref_sequences()
    # sums up second column (read_counts)
    r_clv_ref = sum([s[1] for s in clv_ref_info])

    # 2 - Select and count the reads of the uncleaved sequences in the +ligand round
    unclv_ref_info = db.get_ref_sequences(cleaved_prefix=0)
    # sums up second column (read_counts)
    r_unclv_ref = sum([p[1] for p in unclv_ref_info])

    # 3 - Select and count the reads of the cleaved sequences in the +ligand round
    clv_seq_info = db.get_sequences()  # returns list of tuples with one value?

    for i in range(len(clv_seq_info)):
        # get specific tuple with information of the sequence
        one_seq_info = clv_seq_info[i]
        one_seq = one_seq_info[3]  # get cleaned sequence
        one_ID = one_seq_info[0]  # get key ID

        r_clv = one_seq_info[1]  # get read count

        # find uncleaved reads of the same sequence
        unclv_seq_info = db.get_uncleaved_sequence(cleaned_sequence=one_seq)
        unclv_ID = unclv_seq_info[0]  # get key ID of uncleaved sequence

        r_unclv = unclv_seq_info[1]  # get read count

        # calculate the cleavage fraction
        clvg_frac = (r_clv / r_clv_ref) / \
            ((r_unclv / r_unclv_ref) + (r_clv / r_clv_ref))
