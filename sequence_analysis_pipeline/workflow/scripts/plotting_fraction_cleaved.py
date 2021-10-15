import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from database_interface import DatabaseInterfaceSequences

database_path = "./results/databases/S1_D63_5_database.db"
# database_path = snakemake.input[0]

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"

with DatabaseInterfaceSequences(path=database_path) as db:
    # Retrieve all unique sequences
    id_sequence_rows = db.get(table=TABLE_ID_SEQ, columns=["id"])

    num_ids = len(id_sequence_rows)

    # Id position is a flag for a possible biosensor
    biosensor_array = np.zeros(num_ids+1, dtype=np.int16)

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

fig2, ax2 = plt.subplots(figsize=(14, 8), dpi=125)
ax2.scatter(fraction_negative[non_biosensor_array],
            fraction_positive[non_biosensor_array], s=5, alpha=1, color="#057D54")
ax2.scatter(fraction_negative[biosensor_array],
            fraction_positive[biosensor_array], s=5, alpha=1, color="red")

ax2.plot(np.arange(100), np.arange(100),
         color="#9B0138", label="reference", linestyle="dashed", linewidth=1)


ax2.set_xlim([0, 100])
ax2.set_ylim([0, 100])
ax2.xaxis.set_major_locator(MultipleLocator(10))
ax2.xaxis.set_minor_locator(MultipleLocator(5))
ax2.yaxis.set_major_locator(MultipleLocator(10))
ax2.yaxis.set_minor_locator(MultipleLocator(5))

ax2.set_xlabel("Fraction cleaved (-folate)")
ax2.set_ylabel("Fraction cleaved (+folate)")

# ax2.set_xscale("symlog")
# # ax2.set_yscale("log")
# fig2.show()
plt.show()


# ax2.set_yscale("log")
# ax2.set_ylim((0, 10))
# ax2.spines["top"].set_visible(False)
# ax2.xaxis.set_ticks_position("bottom")

# divider = make_axes_locatable(ax2)
# ax2_lin = divider.append_axes("top", size=2.0, pad=0, sharex=ax2)
# ax2_lin.scatter(fraction_negative[non_biosensor_array],
#                 fraction_positive[non_biosensor_array], s=10, alpha=1)
# ax2_lin.scatter(fraction_negative[biosensor_array],
#                 fraction_positive[biosensor_array], s=5, alpha=1, color="red")
# ax2_lin.set_xscale("linear")
# ax2_lin.set_ylim((10, 90))

# ax2_log2 = divider.append_axes("top", size=2.0, pad=0, sharex=ax2)
# ax2_log2.scatter(fraction_negative[non_biosensor_array],
#                  fraction_positive[non_biosensor_array], s=10, alpha=1)
# ax2_log2.scatter(fraction_negative[biosensor_array],
#                  fraction_positive[biosensor_array], s=5, alpha=1, color="red")
# ax2_log2.set_yscale("log")
# ax2_log2.set_ylim((90, 100))

# # Removes the bottom axis line
# ax2_lin.spines["bottom"].set_visible(False)
# ax2_lin.xaxis.set_ticks_position("top")
# plt.setp(ax2_lin.get_xticklabels(), visible=False)
