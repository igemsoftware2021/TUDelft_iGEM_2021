import numpy as np
import matplotlib.pyplot as plt
from database_interface import DatabaseInterfaceCleanSequences

database_path = "./results/databases/T1_D80_database.db"
# database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    sequences_info = db.query(
        f"SELECT cleaned_sequence, fold_change, fold_change_standard_error FROM {TABLE_NAME} WHERE fold_change IS NOT NULL", fetchall=True)
    uniq_sequences_info = list(set(sequences_info))

print(sequences_info[0])

fold_change_array = np.array(
    [sequence_info[1] for sequence_info in uniq_sequences_info], dtype=np.float64)
fold_change_se_array = np.array(
    [sequence_info[2] for sequence_info in uniq_sequences_info], dtype=np.float64)

min_se = np.amin(fold_change_se_array)

print(min_se)

fold_change_se_array_updated = fold_change_se_array
print(fold_change_se_array_updated)

fig, ax = plt.subplots()
ax.scatter(fold_change_array, fold_change_se_array_updated, s=2, alpha=1)
ax.grid()
# ax.set_ylim([0.0, 0.02])
ax.set_xscale('log')
ax.set_yscale('log')
# fig.savefig("fold_change_vs_fold_change_se.pdf", dpi=300)
plt.show()
fig.show()
