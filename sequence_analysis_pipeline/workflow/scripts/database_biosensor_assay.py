# add scipy.stats to
import scipy.stats
from database_interface import DatabaseInterfaceSequences

# database_path = "./results/databases/S1_D63_15_database.db"
database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

# We need to do one-sided
with DatabaseInterfaceSequences(path=database_path) as db:
    sequences_info = db.query(
        f"SELECT id, sequence, fold_change_estimated_mean, fold_change_standard_error FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)
    uniq_sequences_info = list(set(sequences_info))

count = 0

N = len(uniq_sequences_info)
alpha = 1/N
expected_fold_change = 1.0

sequence_set = set()

with DatabaseInterfaceSequences(path=database_path) as db:
    for sequence_info in uniq_sequences_info:
        rowid = sequence_info[0]
        seq_fold_change_estimated_mean = sequence_info[2]
        seq_fold_change_se = sequence_info[3]
        z_score = (seq_fold_change_estimated_mean -
                   expected_fold_change) / seq_fold_change_se
        p_value = scipy.stats.norm.sf(abs(z_score))  # one-sided
        db.update_column_value(TABLE_NAME, rowid, "p_value", p_value)

        # if seq_fold_change > 1.4:
        #     print(rowid, sequence_info)
        if seq_fold_change_estimated_mean > 3.0 and p_value < 1/N:
            print(rowid, sequence_info)
            db.update_column_value(TABLE_NAME, rowid, "possible_sensor", 1)
            count += 1

            sequence_set.add(
                (sequence_info[1], seq_fold_change_estimated_mean, p_value))
            # if p_value < alpha:
            # count += 1
            # print(sequence_info[0], seq_fold_change, p_value)

for info in sequence_set:
    print(info[0])
