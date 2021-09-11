import numpy as np
from tqdm import tqdm
from database_interface import DatabaseInterfaceSequences
from bootstrapping_v2 import bootstrap_cs_fc_with_replacement


# database_path = snakemake.input[0]
database_path = "results/databases/T1_D80_database_v2.db"

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    info_rows = db.get(table=TABLE_CLEAN_SEQ, columns=[
                       "id", "read_count", "sequence_id", "reference_name", "cleaved_prefix", "ligand_present"])

    # for info_row in info_rows:
    #     if info_row[3] is not None:
    #         print(info_row)

    rowid_to_neg_cs_mean, rowid_to_neg_cs_se, rowid_to_pos_cs_mean, rowid_to_pos_cs_se, rowid_to_fc_mean, rowid_to_fc_se = bootstrap_cs_fc_with_replacement(
        info_rows)

    for rowid in rowid_to_neg_cs_mean.keys():
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="cleavage_fraction_estimated_mean", value=rowid_to_neg_cs_mean[rowid])
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="cleavage_fraction_standard_error", value=rowid_to_neg_cs_se[rowid])

    for rowid in rowid_to_pos_cs_mean.keys():
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="cleavage_fraction_estimated_mean", value=rowid_to_pos_cs_mean[rowid])
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="cleavage_fraction_standard_error", value=rowid_to_pos_cs_se[rowid])

    for rowid in rowid_to_fc_mean.keys():
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="fold_change_estimated_mean", value=rowid_to_fc_mean[rowid])
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="fold_change_standard_error", value=rowid_to_fc_se[rowid])
