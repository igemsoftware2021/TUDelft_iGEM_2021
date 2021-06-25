import sqlite3
from database_interface import DatabaseInterfaceCleanSequences
from tqdm import tqdm

# add right path for the table with the clean sequences
# database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:

    # 1 - Get number of reads of the references with cleaved prefix in +ligand round
    clv_ref_info = db.get_info_ref_sequences(table=TABLE_NAME)
    # sums up second column (read_counts)
    r_clv_ref = sum([s[1] for s in clv_ref_info])

    # 2 - Get number of reads of the references with uncleaved prefix in +ligand round
    unclv_ref_info = db.get_info_ref_sequences(
        table=TABLE_NAME, cleaved_prefix=0)
    # sums up second column (read_counts)
    r_unclv_ref = sum([p[1] for p in unclv_ref_info])

    # 3 - Do the same with the references for the -ligand round
    # cleaved prefix
    clv_ref_info_neg = db.get_info_ref_sequences(
        table=TABLE_NAME, ligand_present=0)
    # sums up second column (read_counts)
    r_clv_ref_neg = sum([k[1] for k in clv_ref_info_neg])

    # uncleaved prefix
    unclv_ref_info_neg = db.get_info_ref_sequences(table=TABLE_NAME,
                                                   cleaved_prefix=0, ligand_present=0)
    # sums up second column (read_counts)
    r_unclv_ref_neg = sum([a[1] for a in unclv_ref_info])

    # 4 - Select and count the reads of the cleaved sequences in the +ligand round
    # returns list of tuples with one value?
    clv_seq_info = db.get_sequences(
        table=TABLE_NAME, cleaved_prefix=1, ligand_present=1)

    for i in tqdm(range(len(clv_seq_info))):
        # get specific tuple with information of the sequence
        one_seq_info = clv_seq_info[i]
        one_seq = one_seq_info[2]  # get cleaned sequence
        one_ID = one_seq_info[0]   # get key ID
        r_clv = one_seq_info[1]    # get read count

        # find uncleaved reads of the same sequence
        unclv_seq_info = db.get_info_sequence(
            table=TABLE_NAME, cleaned_sequence=one_seq, cleaved_prefix=0, ligand_present=1)[0]
        unclv_ID = unclv_seq_info[0]  # get key ID of uncleaved sequence
        r_unclv = unclv_seq_info[1]   # get read count of uncleaved sequence

        # calculate the cleavage fraction in +ligand round
        clvg_frac = (r_clv / r_clv_ref) / \
            ((r_unclv / r_unclv_ref) + (r_clv / r_clv_ref))

        #put in database
        db.update_cleavage_fraction(
            table=TABLE_NAME, rowid=one_ID, cleavage_fraction=clvg_frac)
        db.update_cleavage_fraction(
            table=TABLE_NAME, rowid=unclv_ID, cleavage_fraction=clvg_frac)

        # -ligand round
        # get info of cleaved sequence in negative round
        clv_seq_info_neg = db.get_info_sequence(
            table=TABLE_NAME, cleaned_sequence=one_seq, cleaved_prefix=1, ligand_present=0)[0]
        # get info of uncleaved sequence in negative round
        unclv_seq_info_neg = db.get_info_sequence(
            table=TABLE_NAME, cleaned_sequence=one_seq, cleaved_prefix=0, ligand_present=0)[0]

        # cleaved sequence
        # get ID of cleaved sequence in negative round
        clv_ID_neg = clv_seq_info_neg[0]
        # get read count cleaved sequence in negative round
        r_clv_neg = clv_seq_info_neg[1]

        # uncleaved sequence
        # get ID of uncleaved sequence in negative round
        unclv_ID_neg = unclv_seq_info_neg[0]
        # get read count uncleaved sequence in negative round
        r_unclv_neg = unclv_seq_info_neg[1]

        # calculate the cleavage fraction in -ligand round
        clvg_frac_neg = (r_clv_neg / r_clv_ref_neg) / \
            ((r_unclv_neg / r_unclv_ref_neg) + (r_clv_neg / r_clv_ref_neg))

        #put in database
        db.update_cleavage_fraction(
            table=TABLE_NAME, rowid=clv_ID_neg, cleavage_fraction=clvg_frac_neg)
        db.update_cleavage_fraction(
            table=TABLE_NAME, rowid=unclv_ID_neg, cleavage_fraction=clvg_frac_neg)
