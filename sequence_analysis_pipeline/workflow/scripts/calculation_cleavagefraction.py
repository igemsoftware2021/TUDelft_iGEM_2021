from tqdm import tqdm
from database_interface import DatabaseInterfaceCleanSequences
import calc_helpers

# add right path for the table with the clean sequences
# database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
database_path = snakemake.input[0]
# database_path = "results/databases/T1_D80_database.db"

TABLE_NAME = "clean_sequences"


# TODO first check whether it is even possible for a sequence to have a cleavage fraction


with DatabaseInterfaceCleanSequences(path=database_path) as db:

    # Get all the reference sequence info for certain conditions
    ref_unclvd_neg = db.get_info_ref_sequences(
        TABLE_NAME, cleaved_prefix=0, ligand_present=0)
    ref_clvd_neg = db.get_info_ref_sequences(
        TABLE_NAME, cleaved_prefix=1, ligand_present=0)
    ref_unclvd_pos = db.get_info_ref_sequences(
        TABLE_NAME, cleaved_prefix=0, ligand_present=1)
    ref_clvd_pos = db.get_info_ref_sequences(
        TABLE_NAME, cleaved_prefix=1, ligand_present=1)

    # read count is at second position in every tuple
    r_ref_unclvd_neg = sum(seq_info[1] for seq_info in ref_unclvd_neg)
    r_ref_clvd_neg = sum(seq_info[1] for seq_info in ref_clvd_neg)
    r_ref_unclvd_pos = sum(seq_info[1] for seq_info in ref_unclvd_pos)
    r_ref_clvd_pos = sum(seq_info[1] for seq_info in ref_clvd_pos)

    # # 1 - Get number of reads of the references with cleaved prefix in +ligand round
    # clv_ref_info = db.get_info_ref_sequences(table=TABLE_NAME)
    # # sums up second column (read_counts)
    # r_clv_ref = sum([s[1] for s in clv_ref_info])

    # # 2 - Get number of reads of the references with uncleaved prefix in +ligand round
    # unclv_ref_info = db.get_info_ref_sequences(
    #     table=TABLE_NAME, cleaved_prefix=0)
    # # sums up second column (read_counts)
    # r_unclv_ref = sum([p[1] for p in unclv_ref_info])

    # # 3 - Do the same with the references for the -ligand round
    # # cleaved prefix
    # clv_ref_info_neg = db.get_info_ref_sequences(
    #     table=TABLE_NAME, ligand_present=0)
    # # sums up second column (read_counts)
    # r_clv_ref_neg = sum([k[1] for k in clv_ref_info_neg])

    # # uncleaved prefix
    # unclv_ref_info_neg = db.get_info_ref_sequences(table=TABLE_NAME,
    #                                                cleaved_prefix=0, ligand_present=0)
    # # sums up second column (read_counts)
    # r_unclv_ref_neg = sum([a[1] for a in unclv_ref_info])

    # # 4 - Select and count the reads of the cleaved sequences in the +ligand round
    # # returns list of tuples with one value?
    # clv_seq_info = db.get_sequences(
    #     table=TABLE_NAME, cleaved_prefix=1, ligand_present=1)

    # Fetchall all rows only columns (id, cleaned_sequence, prefix_name, ligand_present)
    rows_info_seq = db.get(table=TABLE_NAME, columns=[
        "id", "cleaned_sequence", "cleaved_prefix", "ligand_present"])

    # Fetchall all rows only columns (id, read_count)
    rows_id_and_read_count = db.get(table=TABLE_NAME, columns=[
        "id", "read_count"])

# Create dictionary where (sequence, cleaved_prefix, ligand_present) is the key
# and the value is the rowid corresponding to this
info_seq_rowid_dict = {}
for row_info in rows_info_seq:
    key_tuple = (row_info[1], row_info[2], row_info[3])
    rowid = row_info[0]
    info_seq_rowid_dict[key_tuple] = rowid

# Dictionary where key is the rowid, and the value is the read count
id_read_count_dict = {row_info[0]: row_info[1]
                      for row_info in rows_id_and_read_count}


with DatabaseInterfaceCleanSequences(path=database_path) as db:

    seq_unclvd_neg = db.query(
        f"SELECT cleaned_sequence, cleaved_prefix, ligand_present FROM {TABLE_NAME} WHERE cleaved_prefix=0 AND ligand_present=0", fetchall=True)

    print("Calculating the cleavage fraction for the sequences where ligand was not present...")
    for seq_unclvd_neg_info in tqdm(seq_unclvd_neg):

        rowid_seq_unclvd_neg = info_seq_rowid_dict.get(seq_unclvd_neg_info)
        rowid_seq_clvd_neg = info_seq_rowid_dict.get(
            (seq_unclvd_neg_info[0], 1, seq_unclvd_neg_info[2]))

        r_seq_unclvd_neg = id_read_count_dict.get(rowid_seq_unclvd_neg)
        r_seq_clvd_neg = id_read_count_dict.get(rowid_seq_clvd_neg)

        if r_seq_unclvd_neg is None or r_seq_clvd_neg is None:
            continue

        # else
        cs_neg = calc_helpers.calc_cleavage_fraction(
            r_seq_clvd_neg, r_ref_clvd_neg, r_seq_unclvd_neg, r_ref_unclvd_neg)

        # Update the cleavage fraction in the database
        db.update_column_value(table=TABLE_NAME, rowid=rowid_seq_unclvd_neg,
                               column_name="cleavage_fraction", value=cs_neg)
        db.update_column_value(table=TABLE_NAME, rowid=rowid_seq_clvd_neg,
                               column_name="cleavage_fraction", value=cs_neg)

    seq_unclvd_pos = db.query(
        f"SELECT cleaned_sequence, cleaved_prefix, ligand_present FROM {TABLE_NAME} WHERE cleaved_prefix=0 AND ligand_present=1", fetchall=True)

    print("Calculating the cleavage fraction for the sequences where ligand was present...")
    for seq_unclvd_pos_info in tqdm(seq_unclvd_pos):

        rowid_seq_unclvd_pos = info_seq_rowid_dict.get(seq_unclvd_pos_info)
        rowid_seq_clvd_pos = info_seq_rowid_dict.get(
            (seq_unclvd_pos_info[0], 1, seq_unclvd_pos_info[2]))

        r_seq_unclvd_pos = id_read_count_dict.get(rowid_seq_unclvd_pos)
        r_seq_clvd_pos = id_read_count_dict.get(rowid_seq_clvd_pos)

        if r_seq_unclvd_pos is None or r_seq_clvd_pos is None:
            continue

        # else
        cs_pos = calc_helpers.calc_cleavage_fraction(
            r_seq_clvd_pos, r_ref_clvd_pos, r_seq_unclvd_pos, r_ref_unclvd_pos)

        # Update the cleavage fraction in the database
        db.update_column_value(table=TABLE_NAME, rowid=rowid_seq_unclvd_pos,
                               column_name="cleavage_fraction", value=cs_pos)
        db.update_column_value(table=TABLE_NAME, rowid=rowid_seq_clvd_pos,
                               column_name="cleavage_fraction", value=cs_pos)
