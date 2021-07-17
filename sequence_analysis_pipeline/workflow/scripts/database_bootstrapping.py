import numpy as np
from tqdm import tqdm
from database_interface import DatabaseInterfaceCleanSequences
import bootstrapping

# database_path = "./results/databases/T1_D80_database.db"
database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

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

    # First retrieve all sequences with a fold change
    sequences_with_fold_change = db.query(
        f"SELECT cleaned_sequence FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)

    # Fetchall all rows only columns (id, cleaned_sequence, prefix_name, ligand_present) where fold change is not NULL
    rows_info_seq = db.query(
        f"SELECT id, cleaned_sequence, cleaved_prefix, ligand_present FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)

    # Fetchall all rows only columns (id, cleavage_fraction) where fold change is not NULL
    rowsid_and_read_count = db.query(
        f"SELECT id, read_count FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)

    # Fetchall all rows where fold change is not NULL only columns (id, cleavage_fraction)
    rowsid_and_k_factor = db.query(
        f"SELECT id, k_factor FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)

# All unique sequences with a fold_change
sequences_with_fold_change = [row_info[0]
                              for row_info in sequences_with_fold_change]

# Create dictionary where (sequence, cleaved_prefix, ligand_present) is the key
# and the value is the rowid corresponding to this
info_seq_rowid_dict = {}
for row_info in rows_info_seq:
    key_tuple = (row_info[1], row_info[2], row_info[3])
    rowid = row_info[0]
    info_seq_rowid_dict[key_tuple] = rowid

# Dictionary where key is the rowid, and the value is the cleavage_fraction
rowid_read_count_dict = {row_info[0]: row_info[1]
                         for row_info in rowsid_and_read_count}

# Dictionary where key is the rowid, and the value is the cleavage_fraction
rowid_k_factor_dict = {row_info[0]: row_info[1]
                       for row_info in rowsid_and_k_factor}

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    for sequence in tqdm(sequences_with_fold_change):

        # Retrieve sequence info condition ligand not present
        rowid_seq_unclvd_neg = info_seq_rowid_dict[(sequence, 0, 0)]
        rowid_seq_clvd_neg = info_seq_rowid_dict[(sequence, 1, 0)]

        # r_ = read_count
        # Retrieve read count of sequence with uncleaved prefix
        r_seq_unclvd_neg = rowid_read_count_dict[rowid_seq_unclvd_neg]
        # Retrieve read count of sequence with cleaved prefix
        r_seq_clvd_neg = rowid_read_count_dict[rowid_seq_clvd_neg]

        # Retrieve sequence info condition ligand present
        rowid_seq_unclvd_pos = info_seq_rowid_dict[(sequence, 0, 1)]
        rowid_seq_clvd_pos = info_seq_rowid_dict[(sequence, 1, 1)]

        # r_ = read_count
        # Retrieve read count of sequence with uncleaved prefix
        r_seq_unclvd_pos = rowid_read_count_dict[rowid_seq_unclvd_pos]
        # Retrieve read count of sequence with cleaved prefix
        r_seq_clvd_pos = rowid_read_count_dict[rowid_seq_clvd_pos]

        # k_factor is the 15th position
        k_factor = rowid_k_factor_dict[rowid_seq_clvd_pos]

        # Now determine the fold change things
        cs_neg_mean, cs_neg_se, cs_pos_mean, cs_pos_se, fold_changes, fold_change_mean, fold_change_se = bootstrapping.bootstrap_fold_change_with_replacement(
            r_seq_clvd_neg, r_seq_unclvd_neg, r_seq_clvd_pos, r_seq_unclvd_pos, r_ref_clvd_neg, r_ref_unclvd_neg, r_ref_clvd_pos, r_ref_unclvd_pos, k=k_factor, num_samples=1000)

        # Update the sequence info in the database for condition ligand not present
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_neg, "cleavage_fraction_estimated_mean", cs_neg_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_neg, "cleavage_fraction_standard_error", cs_neg_se)

        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_neg, "cleavage_fraction_estimated_mean", cs_neg_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_neg, "cleavage_fraction_standard_error", cs_neg_se)

        # Update the sequence info in the database for condition ligand present
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_pos, "cleavage_fraction_estimated_mean", cs_pos_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_pos, "cleavage_fraction_standard_error", cs_pos_se)

        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_pos, "cleavage_fraction_estimated_mean", cs_pos_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_pos, "cleavage_fraction_standard_error", cs_pos_se)

        # First update for ligand condition not present
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_neg, "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_neg, "fold_change_standard_error", fold_change_se)

        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_neg, "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_neg, "fold_change_standard_error", fold_change_se)

        # Then update columns for ligand condition present
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_pos, "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_unclvd_pos, "fold_change_standard_error", fold_change_se)

        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_pos, "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, rowid_seq_clvd_pos, "fold_change_standard_error", fold_change_se)
