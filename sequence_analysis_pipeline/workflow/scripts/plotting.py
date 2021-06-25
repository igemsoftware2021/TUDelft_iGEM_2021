import numpy as np
import matplotlib.pyplot as plt
from database_interface import DatabaseInterfaceCleanSequences

# database_path = "./data/NGS/T1_D80_database.db"
database_path = snakemake.input[0]

TABLE_NAME = "clean_sequences"

with DatabaseInterfaceCleanSequences(path=database_path) as db:
    sequences_info = db.get(
        TABLE_NAME, ["cleaned_sequence", "fold_change", "fold_change_se"])
    uniq_sequences_info = list(set(sequences_info))

fold_change_array = np.array(
    [sequence_info[1] for sequence_info in uniq_sequences_info], dtype=np.float32)
fold_change_se_array = np.array(
    [sequence_info[1] for sequence_info in uniq_sequences_info], dtype=np.float32)

fig, ax = plt.subplots()
ax.plot(fold_change_array, fold_change_se_array)
fig.savefig("fold_change_vs_fold_change_se.pdf", dpi=300)
