import numpy as np
from database_interface import DatabaseInterfaceCleanSequences
import bootstrapping
from workflow.scripts.bootstrapping import bootstrap_cleavage_fraction_with_replacement, bootstrap_fold_change_with_replacement

database_path = "sequence_analysis_pipeline/data/NGS/T1_D80_database.db"
# database_path = snakemake.output[0]

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

        # Ligand not present sequences
        # r_ = read_count
        # Retrieve read count of sequence with uncleaved prefix
        r_seq_unclvd_neg = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=0, ligand_present=0)[0][1]
        # Retrieve read count of sequence with cleaved prefix
        r_seq_clvd_neg = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=1, ligand_present=0)[0][1]

        # Create an array with all zeros
        sample_data_neg = np.zeros(
            r_seq_unclvd_neg+r_seq_clvd_neg, dtype=np.int8)
        # Fill the array with an amount of r_seq_clvd_neg ones.
        sample_data_neg[-r_seq_clvd_neg:] = 1

        cs_neg, cs_neg_sd, cs_neg_5_perc, cs_neg_95_perc = bootstrapping.bootstrap_cleavage_fraction_with_replacement(
            sample_data_neg, r_ref_clvd_neg, r_ref_unclvd_neg)

        # Ligand present sequences
        r_seq_unclvd_pos = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=0, ligand_present=1)[0][1]
        r_seq_clvd_pos = db.get_info_sequence(
            TABLE_NAME, cleaned_sequence=sequence, cleaved_prefix=1, ligand_present=1)[0][1]

        sample_data_pos = np.zeros(
            r_seq_unclvd_pos+r_seq_clvd_pos, dtype=np.int8)
        sample_data_pos[-r_seq_clvd_pos:] = 1

        cs_pos, cs_pos_sd, cs_pos_5_perc, cs_pos_95_perc = bootstrapping.bootstrap_cleavage_fraction_with_replacement(
            sample_data_pos, r_ref_clvd_pos, r_ref_unclvd_pos)

    # TODO change k, read this out from earlier

        # Now determine the fold change things
        fold_change_sd, fold_change_se, fold_change_5_perc, fold_change_95_perc = bootstrapping.bootstrap_fold_change_with_replacement(
            cs_neg, cs_pos, k=1)

    # TODO only store the standard deviations

    # update_column_value(self, table: str, rowid: int, column_name: str, value):
    db.update_column_value(TABLE_NAME, rowid, "fold_change", fold_change_sd)
