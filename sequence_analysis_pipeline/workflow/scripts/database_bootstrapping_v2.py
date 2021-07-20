import numpy as np
from tqdm import tqdm
from database_interface import DatabaseInterfaceCleanSequences
import calc_helpers


def create_sample_and_probabilities(id_read_count: dict):
    total = 0
    for key, value in id_read_count.items():
        total += value

    # The ids in the database start at 1
    sample = np.arange(1, len(dict.keys()), dtype=np.int32)
    probabilities = np.zeros(len(dict.keys()), dtype=np.float64)

    i = 0
    for key, value in id_read_count.items():
        if key == i+1:
            probabilities[i] = value/total
        else:
            raise ValueError(
                "Check ids in the database, something is going wrong. key should be equal to i+1")

    return sample, probabilities


# Number of bootstrap samples
NUM_SAMPLES = 1000

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

# cleavage fractions dict
cleavage_fractions_dict = {}

# Fold changes dict
fold_changes_dict = {}

# Create dictionary where (sequence, cleaved_prefix, ligand_present) is the key
# and the value is the rowid corresponding to this
info_seq_rowid_dict = {}
for row_info in rows_info_seq:
    key_tuple = (row_info[1], row_info[2], row_info[3])
    rowid = row_info[0]
    info_seq_rowid_dict[key_tuple] = rowid

    # Fill the cleavage_fraction and fold changes dict
    cleavage_fractions_dict[(row_info[0], row_info[2])] = []
    fold_changes_dict[row_info[0]] = []

# Dictionary where key is the rowid, and the value is the cleavage_fraction
rowid_read_count_dict = {row_info[0]: row_info[1]
                         for row_info in rowsid_and_read_count}

total_read_counts = 0
for key, value in rowid_read_count_dict.items():
    total_read_counts += value

# Dictionary where key is the rowid, and the value is the cleavage_fraction
rowid_k_factor_dict = {row_info[0]: row_info[1]
                       for row_info in rowsid_and_k_factor}

# Create the sample array and the probabilities array
sample, probabilities = create_sample_and_probabilities(rowid_read_count_dict)


for i in tqdm(range(NUM_SAMPLES)):
    bootstrap_sample = np.random.choice(
        sample, size=total_read_counts, dtype=np.int32)

    unique, counts = np.unique(bootstrap_sample, return_counts=True)
    seq_counts = dict(zip(unique, counts))

    for sequence in sequences_with_fold_change:
        # Retrieve rowids for sequence condition ligand not present
        rowid_seq_unclvd_neg = info_seq_rowid_dict[(sequence, 0, 0)]
        rowid_seq_clvd_neg = info_seq_rowid_dict[(sequence, 1, 0)]

        # Retrieve rowids for sequence condition ligand present
        rowid_seq_unclvd_pos = info_seq_rowid_dict[(sequence, 0, 1)]
        rowid_seq_clvd_pos = info_seq_rowid_dict[(sequence, 1, 1)]

        # Calculate the read counts
        sample_r_unclvd_s_neg = seq_counts[rowid_seq_unclvd_neg]
        sample_r_clvd_s_neg = seq_counts[rowid_seq_clvd_neg]

        sample_r_unclvd_s_pos = seq_counts[rowid_seq_unclvd_pos]
        sample_r_clvd_s_pos = seq_counts[rowid_seq_clvd_pos]

        cs_neg = calc_helpers.calc_cleavage_fraction(
            sample_r_clvd_s_neg, r_ref_clvd_neg, sample_r_unclvd_s_neg, r_ref_unclvd_neg)

        cs_pos = calc_helpers.calc_cleavage_fraction(
            sample_r_clvd_s_pos, r_ref_clvd_pos, sample_r_unclvd_s_pos, r_ref_unclvd_pos)

        fold_change = calc_helpers.calc_fold_change(
            cs_pos, cs_neg, k=rowid_k_factor_dict[rowid_seq_unclvd_neg])

        # Fill the cleavage fractions dict
        temp_list = cleavage_fractions_dict.get((sequence, 0), [])
        temp_list.append(cs_neg)
        cleavage_fractions_dict[(sequence, 0)] = temp_list
        temp_list = cleavage_fractions_dict.get((sequence, 1), [])
        temp_list.append(cs_pos)
        cleavage_fractions_dict[(sequence, 1)] = temp_list

        # Fill the fold changes dict
        temp_list = fold_changes_dict.get(sequence, [])
        temp_list.append(fold_change)
        fold_changes_dict[sequence] = temp_list

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    for sequence in tqdm(sequences_with_fold_change):

        # Retrieve rowid sequence info condition ligand not present
        rowid_seq_unclvd_neg = info_seq_rowid_dict[(sequence, 0, 0)]
        rowid_seq_clvd_neg = info_seq_rowid_dict[(sequence, 1, 0)]

        # Retrieve rowid sequence info condition ligand present
        rowid_seq_unclvd_pos = info_seq_rowid_dict[(sequence, 0, 1)]
        rowid_seq_clvd_pos = info_seq_rowid_dict[(sequence, 1, 1)]

        cs_neg = cleavage_fractions_dict[(sequence, 0)]

        # Calculate estimated mean
        cs_neg_mean = np.mean(cs_neg)

        # Calculate standard deviation of the bootstrapped cleavage fractions
        cs_neg_se = calc_helpers.calc_sample_standard_deviation(cs_neg)

        cs_pos = cleavage_fractions_dict[(sequence, 1)]

        # Calculate estimated mean
        cs_pos_mean = np.mean(cs_pos)

        # Calculate standard deviation of the bootstrapped cleavage fractions
        cs_pos_se = calc_helpers.calc_sample_standard_deviation(cs_pos)

        fold_changes = fold_changes_dict[sequence]

        fold_change_mean = np.mean(fold_changes)

        fold_change_se = calc_helpers.calc_sample_standard_deviation(
            fold_changes)

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
