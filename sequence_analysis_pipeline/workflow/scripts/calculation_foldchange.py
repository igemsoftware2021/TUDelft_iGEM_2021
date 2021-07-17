import numpy as np
from tqdm import tqdm
from database_interface import DatabaseInterfaceCleanSequences
import calc_helpers


# add right path for the table with the clean sequences
# database_path = "sequence_analysis_pipeline/data/NGS/processed/S1_D80_database.db"
database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    # Fetchall all rows only columns (id, cleaned_sequence, prefix_name, ligand_present)
    rows_info_seq = db.get(table=TABLE_NAME, columns=[
        "id", "cleaned_sequence", "cleaved_prefix", "ligand_present"])

    # Fetchall all rows only columns (id, cleavage_fraction)
    rows_id_and_cs = db.get(table=TABLE_NAME, columns=[
        "id", "cleavage_fraction"])

# Create dictionary where (sequence, cleaved_prefix, ligand_present) is the key
# and the value is the rowid corresponding to this
info_seq_rowid_dict = {}
for row_info in rows_info_seq:
    key_tuple = (row_info[1], row_info[2], row_info[3])
    rowid = row_info[0]
    info_seq_rowid_dict[key_tuple] = rowid

# Dictionary where key is the rowid, and the value is the cleavage_fraction
id_cs_dict = {row_info[0]: row_info[1]
              for row_info in rows_id_and_cs}


# 1 - calculate the fold-change for k=1

# temporary storage for the calculated fold-change with k=1
fc_list = []

with DatabaseInterfaceCleanSequences(path=database_path) as db:

    seq_unclvd = db.query(
        f"SELECT cleaned_sequence, cleaved_prefix FROM {TABLE_NAME} WHERE cleaved_prefix=0", fetchall=True)

    for seq_info in tqdm(seq_unclvd):

        rowid_cs_neg_unclvd = info_seq_rowid_dict.get(
            (seq_info[0], seq_info[1], 0))
        rowid_cs_pos_unclvd = info_seq_rowid_dict.get(
            (seq_info[0], seq_info[1], 1))

        cs_neg = id_cs_dict.get(rowid_cs_neg_unclvd)
        cs_pos = id_cs_dict.get(rowid_cs_pos_unclvd)

        if cs_neg is None or cs_pos is None:
            continue

        # calculate fold-change with k=1
        fc_temp = calc_helpers.calc_fold_change(cs_pos, cs_neg, k=1.0)
        fc_list.append(fc_temp)

# 2 - Find median value
# sort the list with temporary fold change values
median_k = np.median(fc_list, overwrite_input=True)

# 3 - Calculate new k with median fold-change value
k_new = 1 / median_k
with DatabaseInterfaceCleanSequences(path=database_path) as db:

    db.query(f"UPDATE {TABLE_NAME} SET k_factor={k_new}")

    seq_unclvd = db.query(
        f"SELECT cleaned_sequence, cleaved_prefix FROM {TABLE_NAME} WHERE cleaved_prefix=0", fetchall=True)

    for seq_info in tqdm(seq_unclvd):

        rowid_cs_neg_unclvd = info_seq_rowid_dict.get(
            (seq_info[0], seq_info[1], 0))
        rowid_cs_pos_unclvd = info_seq_rowid_dict.get(
            (seq_info[0], seq_info[1], 1))

        rowid_cs_neg_clvd = info_seq_rowid_dict.get(
            (seq_info[0], 1, 0))
        rowid_cs_pos_clvd = info_seq_rowid_dict.get(
            (seq_info[0], 1, 1))

        cs_neg = id_cs_dict.get(rowid_cs_neg_unclvd)
        cs_pos = id_cs_dict.get(rowid_cs_pos_unclvd)

        if cs_neg is None or cs_pos is None:
            continue

        # calculate fold-change with k=k_new
        fc = calc_helpers.calc_fold_change(cs_pos, cs_neg, k=k_new)

        db.update_column_value(table=TABLE_NAME, rowid=rowid_cs_neg_unclvd,
                               column_name="fold_change", value=fc)
        db.update_column_value(table=TABLE_NAME, rowid=rowid_cs_neg_clvd,
                               column_name="fold_change", value=fc)
        db.update_column_value(table=TABLE_NAME, rowid=rowid_cs_pos_unclvd,
                               column_name="fold_change", value=fc)
        db.update_column_value(table=TABLE_NAME, rowid=rowid_cs_pos_clvd,
                               column_name="fold_change", value=fc)

    # for i in tqdm(range(len(sequences_data))):
    #     seq = sequences_data[i][2]            # get the sequence
    #     cs_pos = sequences_data[i][9]         # get cleavage fraction
    #     seq_ID = sequences_data[i][0]

    #     # get info of the same sequence but then in the ligand NOT present round, with cleaved_prefix, cause for biosensors, cleaved is definitaly present in -round
    #     sequences_data_neg = db.get_info_sequence(
    #         table=TABLE_NAME, cleaned_sequence=seq, cleaved_prefix=1, ligand_present=0)
    #     cs_neg = sequences_data_neg[0][9]     # get cleavage fraction
    #     # get sequence ID of negative round
    #     seq_ID_neg = sequences_data_neg[0][0]

    #     # calculate fold-change with k_new
    #     fc = k_new * ((1-cs_pos) / (1-cs_neg))

    #     # store the fold change in the database
    #     db.update_fold_change(table=TABLE_NAME, rowid=seq_ID, fold_change=fc)
    #     db.update_fold_change(
    #         table=TABLE_NAME, rowid=seq_ID_neg, fold_change=fc)
