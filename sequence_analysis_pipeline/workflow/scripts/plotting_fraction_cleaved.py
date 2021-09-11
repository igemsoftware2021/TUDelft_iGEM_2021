import numpy as np
import matplotlib.pyplot as plt
from database_interface import DatabaseInterfaceSequences

database_path = "./results/databases/T1_D80_database_v2.db"
# database_path = snakemake.input[0]

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    # Retrieve all unique sequences
    id_sequence_rows = db.get(table=TABLE_ID_SEQ, columns=["id"])

    num_ids = len(id_sequence_rows)

    biosensor_array = np.zeros(num_ids+1, dtype=np.int8)

    # First do negative ligand
    fraction_negative = np.zeros(num_ids+1, dtype=np.float64)
    for id_sequence_row in id_sequence_rows:
        id_sequence = id_sequence_row[0]
        info_cleaved = db.retrieve_info_sequence_id(
            table=TABLE_CLEAN_SEQ, sequence_id=id_sequence, cleaved_prefix=1, ligand_present=0)
        if len(info_cleaved) == 0:
            continue
        else:
            cleavage_fraction = info_cleaved[0][10]
            fraction_negative[id_sequence] = cleavage_fraction

            biosensor = info_cleaved[0][12]
            if biosensor == 1:
                biosensor_array[id_sequence] = 1
                print(id_sequence)

    # Then do for ligand present
    fraction_positive = np.zeros(num_ids+1, dtype=np.float64)
    for id_sequence_row in id_sequence_rows:
        id_sequence = id_sequence_row[0]
        info_cleaved = db.retrieve_info_sequence_id(
            table=TABLE_CLEAN_SEQ, sequence_id=id_sequence, cleaved_prefix=1, ligand_present=1)
        if len(info_cleaved) == 0:
            continue
        else:
            cleavage_fraction = info_cleaved[0][10]
            fraction_positive[id_sequence] = cleavage_fraction

fraction_negative = fraction_negative * 100
fraction_positive = fraction_positive * 100

non_biosensor_array = biosensor_array == 0
biosensor_array = biosensor_array == 1

fig2, ax2 = plt.subplots()
ax2.scatter(fraction_negative[non_biosensor_array],
            fraction_positive[non_biosensor_array], s=2, alpha=1)
ax2.scatter(fraction_negative[biosensor_array],
            fraction_positive[biosensor_array], s=5, alpha=1, color="red")
ax2.set_xlim([0, 100])
ax2.set_ylim([0, 100])
# ax2.set_xscale("log")
# ax2.set_yscale("log")
fig2.show()
plt.show()
