from database_interface import DatabaseInterfaceSequences
from calc_helpers import calc_cs_and_fc_metrics

# database_path = snakemake.input[0]
database_path = "results/databases/T1_D80_database_v2.db"

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    info_rows = db.get(table=TABLE_CLEAN_SEQ, columns=[
                       "id", "read_count", "sequence_id", "reference_name", "cleaved_prefix", "ligand_present"])
    rowid_to_neg_cs, rowid_to_pos_cs, rowid_to_fc = calc_cs_and_fc_metrics(
        info_rows)

    # Store the negative cleavage fractions
    for rowid in rowid_to_neg_cs.keys():
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="cleavage_fraction", value=rowid_to_neg_cs[rowid])
        info_row = db.get_info_rowid(table=TABLE_CLEAN_SEQ, rowid=rowid)
        seq_id = info_row[0][3]
        db.update_column_value(
            table=TABLE_ID_SEQ, rowid=seq_id, column_name="negative_cs", value=rowid_to_neg_cs[rowid])

    # Store the positive cleavage fractions
    for rowid in rowid_to_pos_cs.keys():
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="cleavage_fraction", value=rowid_to_pos_cs[rowid])
        info_row = db.get_info_rowid(table=TABLE_CLEAN_SEQ, rowid=rowid)
        seq_id = info_row[0][3]
        db.update_column_value(
            table=TABLE_ID_SEQ, rowid=seq_id, column_name="positive_cs", value=rowid_to_pos_cs[rowid])

    # Store the fold change values
    for rowid in rowid_to_fc.keys():
        db.update_column_value(table=TABLE_CLEAN_SEQ, rowid=rowid,
                               column_name="fold_change", value=rowid_to_fc[rowid])
        info_row = db.get_info_rowid(table=TABLE_CLEAN_SEQ, rowid=rowid)
        seq_id = info_row[0][3]
        db.update_column_value(
            table=TABLE_ID_SEQ, rowid=seq_id, column_name="fold_change", value=rowid_to_fc[rowid])
