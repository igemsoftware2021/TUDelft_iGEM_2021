import numpy as np
from database_interface import DatabaseInterfaceCleanSequences
import bootstrapping

# database_path = "./data/NGS/T1_D80_database.db"
database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:

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

    # First retrieve all sequences
    sequences_all = db.get(TABLE_NAME, columns=["cleaned_sequence"])
    # Remove duplicate sequences
    sequences_unique = list(set(sequences_all))

    for sequence in sequences_unique:

        # Retrieve sequence info condition ligand not present
        seq_unclvd_neg = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=0, ligand_present=0)[0]
        seq_clvd_neg = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=1, ligand_present=0)[0]

        # r_ = read_count
        # Retrieve read count of sequence with uncleaved prefix
        r_seq_unclvd_neg = seq_unclvd_neg[1]
        # Retrieve read count of sequence with cleaved prefix
        r_seq_clvd_neg = seq_unclvd_neg[1]

        # Create an array with all zeros
        sample_data_neg = np.zeros(
            r_seq_unclvd_neg+r_seq_clvd_neg, dtype=np.int8)
        # Fill the array with an amount of r_seq_clvd_neg ones.
        sample_data_neg[-r_seq_clvd_neg:] = 1

        cs_neg, cs_neg_mean, cs_neg_sd = bootstrapping.bootstrap_cleavage_fraction_with_replacement(
            sample_data_neg, r_ref_clvd_neg, r_ref_unclvd_neg)

        # Update the sequence info in the database for condition ligand not present
        db.update_column_value(
            TABLE_NAME, seq_unclvd_neg[0], "cleavage_fraction_estimated_mean", cs_neg_mean)
        db.update_column_value(
            TABLE_NAME, seq_unclvd_neg[0], "cleavage_fraction_standard_deviation", cs_neg_sd)

        db.update_column_value(
            TABLE_NAME, seq_clvd_neg[0], "cleavage_fraction_estimated_mean", cs_neg_mean)
        db.update_column_value(
            TABLE_NAME, seq_clvd_neg[0], "cleavage_fraction_standard_deviation", cs_neg_sd)

        # Retrieve sequence info condition ligand present
        seq_unclvd_pos = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=0, ligand_present=1)[0]
        seq_clvd_pos = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=1, ligand_present=1)[0]

        # Retrieve read count of sequence with uncleaved prefix
        r_seq_unclvd_pos = seq_unclvd_pos[1]
        # Retrieve read count of sequence with cleaved prefix
        r_seq_clvd_pos = seq_unclvd_pos[1]

        sample_data_pos = np.zeros(
            r_seq_unclvd_pos+r_seq_clvd_pos, dtype=np.int8)
        sample_data_pos[-r_seq_clvd_pos:] = 1

        cs_pos, cs_pos_mean, cs_pos_sd = bootstrapping.bootstrap_cleavage_fraction_with_replacement(
            sample_data_pos, r_ref_clvd_pos, r_ref_unclvd_pos)

        # Update the sequence info in the database for condition ligand present
        db.update_column_value(
            TABLE_NAME, seq_unclvd_pos[0], "cleavage_fraction_estimated_mean", cs_pos_mean)
        db.update_column_value(
            TABLE_NAME, seq_unclvd_pos[0], "cleavage_fraction_standard_deviation", cs_pos_sd)

        db.update_column_value(
            TABLE_NAME, seq_clvd_pos[0], "cleavage_fraction_estimated_mean", cs_pos_mean)
        db.update_column_value(
            TABLE_NAME, seq_clvd_pos[0], "cleavage_fraction_standard_deviation", cs_pos_sd)

        # k_factor is the 15th position
        k_factor = seq_clvd_pos[14]

        # Now determine the fold change things
        fold_change_mean, fold_change_sd, fold_change_se = bootstrapping.bootstrap_fold_change_with_replacement(
            cs_neg, cs_pos, k=k_factor)

        # First update for ligand condition not present
        db.update_column_value(
            TABLE_NAME, seq_unclvd_neg[0], "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, seq_unclvd_neg[0], "fold_change_standard_deviation", fold_change_sd)
        db.update_column_value(
            TABLE_NAME, seq_unclvd_neg[0], "fold_change_standard_error", fold_change_se)

        db.update_column_value(
            TABLE_NAME, seq_clvd_neg[0], "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, seq_clvd_neg[0], "fold_change_standard_deviation", fold_change_sd)
        db.update_column_value(
            TABLE_NAME, seq_clvd_neg[0], "fold_change_standard_error", fold_change_se)

        # Then update columns for ligand condition present
        db.update_column_value(
            TABLE_NAME, seq_unclvd_pos[0], "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, seq_unclvd_pos[0], "fold_change_standard_deviation", fold_change_sd)
        db.update_column_value(
            TABLE_NAME, seq_unclvd_pos[0], "fold_change_standard_error", fold_change_se)

        db.update_column_value(
            TABLE_NAME, seq_clvd_pos[0], "fold_change_estimated_mean", fold_change_mean)
        db.update_column_value(
            TABLE_NAME, seq_clvd_pos[0], "fold_change_standard_deviation", fold_change_sd)
        db.update_column_value(
            TABLE_NAME, seq_clvd_pos[0], "fold_change_standard_error", fold_change_se)
