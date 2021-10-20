import numpy as np
from tqdm import tqdm
import traceback
from calc_helpers import calc_cs_and_fc_metrics, calc_sample_standard_deviation
np.seterr(divide="raise")


def bootstrap_cs_fc_with_replacement(info_rows, num_samples=1000, sample_size=None):
    # info rows looks as follows [(rowid, read_count, sequence_id, reference_name, cleaved_prefix, ligand_present), ..., ()]

    # First create rowid array, where at every index is a rowid integer
    rc_total = sum(info_row[1] for info_row in info_rows)

    rowid_array = np.zeros(rc_total, dtype=np.int32)

    start_idx = 0
    for info_row in info_rows:
        rc = info_row[1]
        rowid_array[start_idx:start_idx+rc] = info_row[0]
        start_idx += rc

    # The dictionaries to store all the information in
    rowid_to_neg_cs_array = dict()
    rowid_to_pos_cs_array = dict()
    rowid_to_fc_array = dict()  # fc = fold change

    for _ in tqdm(range(num_samples)):

        flag = True
        while flag:
            try:
                if sample_size is None:
                    random_idx = np.random.randint(
                        1, high=rowid_array.shape[0], size=rc_total)
                else:
                    random_idx = np.random.randint(
                        1, high=rowid_array.shape[0], size=sample_size)

                random_rowid = rowid_array[random_idx]

                unique, counts = np.unique(random_rowid, return_counts=True)
                rowid_to_rc = dict(zip(unique, counts))

                updated_info_rows = []
                for info_row in info_rows:
                    rowid = info_row[0]
                    rc_new = rowid_to_rc.get(rowid, 0)
                    updated_info_row = (
                        info_row[0], rc_new, info_row[2], info_row[3], info_row[4], info_row[5])
                    updated_info_rows.append(updated_info_row)

                # Calculate the cleavage fractions and fold change
                rowid_to_neg_cs, rowid_to_pos_cs, rowid_to_fc = calc_cs_and_fc_metrics(
                    updated_info_rows)
                flag = False
            except (RuntimeWarning, FloatingPointError) as e:
                # traceback.print_exc()
                # print(e)
                flag = True
                print("again")

        # Append the calculated value to the array in the dictionary
        for rowid in rowid_to_neg_cs.keys():
            value = rowid_to_neg_cs_array.get(rowid, [])
            value.append(rowid_to_neg_cs[rowid])
            rowid_to_neg_cs_array[rowid] = value

        # Append the calculated value to the array in the dictionary
        for rowid in rowid_to_pos_cs.keys():
            value = rowid_to_pos_cs_array.get(rowid, [])
            value.append(rowid_to_pos_cs[rowid])
            rowid_to_pos_cs_array[rowid] = value

        # Append the calculated value to the array in the dictionary
        for rowid in rowid_to_fc.keys():
            value = rowid_to_fc_array.get(rowid, [])
            value.append(rowid_to_fc[rowid])
            rowid_to_fc_array[rowid] = value

    # The dictionaries to store all the information in which will be returned
    rowid_to_neg_cs_mean = dict()
    rowid_to_neg_cs_se = dict()

    for rowid in rowid_to_neg_cs_array.keys():
        neg_cs_array = np.array(rowid_to_neg_cs_array[rowid], dtype=np.float64)
        rowid_to_neg_cs_mean[rowid] = np.mean(neg_cs_array)
        rowid_to_neg_cs_se[rowid] = calc_sample_standard_deviation(
            neg_cs_array)

    # The dictionaries to store all the information in which will be returned
    rowid_to_pos_cs_mean = dict()
    rowid_to_pos_cs_se = dict()

    for rowid in rowid_to_pos_cs_array.keys():
        pos_cs_array = np.array(rowid_to_pos_cs_array[rowid], dtype=np.float64)
        rowid_to_pos_cs_mean[rowid] = np.mean(pos_cs_array)
        rowid_to_pos_cs_se[rowid] = calc_sample_standard_deviation(
            pos_cs_array)

    # The dictionaries to store all the information in which will be returned
    rowid_to_fc_mean = dict()
    rowid_to_fc_se = dict()

    for rowid in rowid_to_fc_array.keys():
        fc_array = np.array(rowid_to_fc_array[rowid], dtype=np.float64)
        rowid_to_fc_mean[rowid] = np.mean(fc_array)
        rowid_to_fc_se[rowid] = calc_sample_standard_deviation(
            fc_array)

    return rowid_to_neg_cs_mean, rowid_to_neg_cs_se, rowid_to_pos_cs_mean, rowid_to_pos_cs_se, rowid_to_fc_mean, rowid_to_fc_se
