import statistics

from workflow.scripts.calc_helpers import calc_cleavage_fraction, calc_fold_change


def calculate_cs_and_fc_metrics(info_rows):
    """Function calculates cleavage fraction and fold change and returns three dictionaries:
    - id to negative cleavage fraction
    - id to positive cleavage fraction
    - id to fold change
    """
    # info rows looks as follows [(rowid, read_count, sequence_id, reference_name, cleaved_prefix, ligand_present), ..., ()]

    # The dictionaries to store all the information in
    rowid_to_neg_cs = dict()
    rowid_to_pos_cs = dict()
    rowid_to_fc = dict()  # fc = fold change

    # Store all the fold change values to change them later after calculating the median
    fc_store = []

    # Store all the sequence ids in a set
    seq_ids = set(info_row[2] for info_row in info_rows)

    # Reference sequence ids
    ref_seq_ids = set(info_row[2]
                      for info_row in info_rows if info_row[3] is not None)

    # Create an dictionary to access all the data you need easily
    info_to_rowid = dict()
    info_to_read_count = dict()
    for info_row in info_rows:
        # (sequence_id, cleaved_prefix, ligand_present) = rowid
        info_to_rowid[(info_row[2], info_row[4], info_row[5])] = info_row[0]
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
            (ref_seq_id, 1, 1), default=0)
        rc_ref_seq_unclvd_pos = info_to_read_count.get(
            (ref_seq_id, 0, 1), default=0)
        rc_ref_seq_clvd_neg = info_to_read_count.get(
            (ref_seq_id, 1, 0), default=0)
        rc_ref_seq_unclvd_neg = info_to_read_count.get(
            (ref_seq_id, 0, 0), default=0)

    # loop over all the sequences and calculate their respective neg cs, pos cs and fc
    for seq_id in seq_ids:

        # First calculate the negative cleavage fraction
        rc_seq_clvd_neg = info_to_read_count.get((seq_id, 1, 0), default=None)
        rc_seq_unclvd_neg = info_to_read_count.get(
            (seq_id, 0, 0), default=None)

        neg_cs_exists = False
        if rc_seq_clvd_neg is not None and rc_seq_unclvd_neg is not None:
            neg_cs = calc_cleavage_fraction(
                rc_seq_clvd_neg, rc_ref_seq_clvd_neg, rc_seq_unclvd_neg, rc_ref_seq_unclvd_neg)
            rowid_seq_clvd_neg = info_to_rowid.get(
                (seq_id, 1, 0), default=None)
            rowid_seq_unclvd_neg = info_to_rowid.get(
                (seq_id, 0, 0), default=None)
            rowid_to_neg_cs[rowid_seq_clvd_neg] = neg_cs
            rowid_to_neg_cs[rowid_seq_unclvd_neg] = neg_cs
            neg_cs_exists = True

        # Then calculate the positive cleavage fraction
        rc_seq_clvd_pos = info_to_read_count.get((seq_id, 1, 1), default=None)
        rc_seq_unclvd_pos = info_to_read_count.get(
            (seq_id, 0, 1), default=None)

        pos_cs_exists = False
        if rc_seq_clvd_pos is not None and rc_seq_unclvd_pos is not None:
            pos_cs = calc_cleavage_fraction(
                rc_seq_clvd_pos, rc_ref_seq_clvd_pos, rc_seq_unclvd_pos, rc_ref_seq_unclvd_pos)
            rowid_seq_clvd_pos = info_to_rowid.get(
                (seq_id, 1, 0), default=None)
            rowid_seq_unclvd_pos = info_to_rowid.get(
                (seq_id, 0, 0), default=None)
            rowid_to_pos_cs[rowid_seq_clvd_pos] = pos_cs
            rowid_to_pos_cs[rowid_seq_unclvd_pos] = pos_cs
            pos_cs_exists = True

        # Finally calculate the fold change if possible
        if neg_cs_exists and pos_cs_exists:
            fc = calc_fold_change(cs_pos=pos_cs, cs_neg=neg_cs, k=1.0)
            fc_store.append(fc)
            rowid_to_fc[rowid_seq_clvd_neg] = fc
            rowid_to_fc[rowid_seq_unclvd_neg] = fc
            rowid_to_fc[rowid_seq_clvd_pos] = fc
            rowid_to_fc[rowid_seq_unclvd_pos] = fc

    # Determine k-factor
    fc_median = statistics.median(fc_store)
    k_new = 1 / fc_median

    # Recalculate the new fold changes
    for seq_id in rowid_to_fc.keys():
        rowid_to_fc[seq_id] = rowid_to_fc[seq_id] * k_new

    return rowid_to_neg_cs, rowid_to_pos_cs, rowid_to_fc
