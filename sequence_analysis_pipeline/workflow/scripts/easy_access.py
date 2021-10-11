from database_interface import DatabaseInterfaceSequences
from collections import Counter
database_path = "./results/databases/S1_D63_database.db"
# database_path = snakemake.input[0]

TABLE_ID_SEQ = "clean_sequences"
# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    # results = db.query(
    #     f"SELECT id, read_count, sequence_id, sequence, reference_name, cleaved_prefix, ligand_present FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)
    # id_sequence_rows = db.get(
    #     table=TABLE_ID_SEQ, columns=["sequence_id"])

    # seq_ids = [seq[0] for seq in id_sequence_rows]
    # print(Counter(seq_ids))

    # print([seq for seq in id_sequence_rows if seq[1] != "Z1"])

    # num_ids = len(id_sequence_rows)
    # for rowid in range(num_ids):
    #     db.update_column_value(TABLE_ID_SEQ, rowid, "possible_sensor", 0)

    # for result in results:
    #     print(result)

    info_rows = db.get(table=TABLE_CLEAN_SEQ, columns=[
                       "id", "read_count", "sequence_id", "reference_name", "cleaved_prefix", "ligand_present"])
    # Reference sequence ids
    ref_seq_ids = set(info_row[2]
                      for info_row in info_rows if info_row[3] is not None)

    info_to_read_count = dict()
    for info_row in info_rows:
        # (sequence_id, cleaved_prefix, ligand_present) = read_count
        info_to_read_count[(info_row[2], info_row[4],
                            info_row[5])] = info_row[1]

    # rc = read count
    # pos = ligand present, neg = ligand not present
    rc_ref_seq_clvd_pos = 0
    rc_ref_seq_unclvd_pos = 0
    rc_ref_seq_clvd_neg = 0
    rc_ref_seq_unclvd_neg = 0

    for ref_seq_id in ref_seq_ids:
        rc_ref_seq_clvd_pos += info_to_read_count.get(
            (ref_seq_id, 1, 1), 0)
        rc_ref_seq_unclvd_pos += info_to_read_count.get(
            (ref_seq_id, 0, 1), 0)
        rc_ref_seq_clvd_neg += info_to_read_count.get(
            (ref_seq_id, 1, 0), 0)
        rc_ref_seq_unclvd_neg += info_to_read_count.get(
            (ref_seq_id, 0, 0), 0)

    print(f"clvd_pos: {rc_ref_seq_clvd_pos}, unclvd_pos: {rc_ref_seq_unclvd_pos}, clvd_neg: {rc_ref_seq_clvd_neg}, unclvd_neg: {rc_ref_seq_unclvd_neg}")
