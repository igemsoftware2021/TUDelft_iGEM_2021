# add scipy.stats to
import scipy.stats
from database_interface import DatabaseInterfaceCleanSequences

database_path = "./results/databases/T1_D80_database.db"
# database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

# We need to do one-sided
with DatabaseInterfaceCleanSequences(path=database_path) as db:
    sequences_info = db.query(
        f"SELECT cleaned_sequence, fold_change, fold_change_standard_error FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)
    uniq_sequences_info = list(set(sequences_info))

count = 0

N = len(uniq_sequences_info)
alpha = 1/N
expected_fold_change = 1.0

for sequence_info in uniq_sequences_info:
    seq_fold_change = sequence_info[1]
    if seq_fold_change > 3.0:
        seq_fold_change_se = sequence_info[2]
        z_score = (seq_fold_change - expected_fold_change) / seq_fold_change_se
        p_value = scipy.stats.norm.sf(abs(z_score))  # one-sided
        count += 1
        print(sequence_info[0], seq_fold_change, p_value)
        # if p_value < alpha:
        # count += 1
        # print(sequence_info[0], seq_fold_change, p_value)
print(count)
