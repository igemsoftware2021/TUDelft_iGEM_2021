import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from database_interface import DatabaseInterfaceSequences


# database_path = "./results/databases/S2_D63_database.db"
database_path = snakemake.input[0]

# The table that links an integer to a sequence
TABLE_ID_SEQ = "id_sequence"

# The table with the added read counts of sequences
TABLE_CLEAN_SEQ = "clean_sequences"


with DatabaseInterfaceSequences(path=database_path) as db:
    # sequences_info = db.query(
    #     f"SELECT sequence, fold_change, fold_change_standard_error, p_value FROM {TABLE_CLEAN_SEQ} WHERE fold_change IS NOT NULL", fetchall=True)
    # uniq_sequences_info = list(set(sequences_info))

    # Retrieve all unique sequences
    id_sequence_rows = db.get(table=TABLE_ID_SEQ, columns=["id"])

    num_ids = len(id_sequence_rows)

    biosensor_array = np.zeros(num_ids+1, dtype=np.int8)
    fc_array = np.zeros(num_ids+1, dtype=np.float64)
    fc_se_array = np.zeros(num_ids+1, dtype=np.float64)
    p_value_array = np.zeros(num_ids+1, dtype=np.float64)

    # First do negative ligand
    for id_sequence_row in id_sequence_rows:
        id_sequence = id_sequence_row[0]
        info_cleaved = db.retrieve_info_sequence_id(
            table=TABLE_CLEAN_SEQ, sequence_id=id_sequence, cleaved_prefix=1, ligand_present=0)
        if len(info_cleaved) == 0:
            continue
        else:
            fc = info_cleaved[0][11]
            fc_se = info_cleaved[0][17]
            biosensor = info_cleaved[0][12]
            p_value = info_cleaved[0][18]

            fc_array[id_sequence] = fc
            fc_se_array[id_sequence] = fc_se
            p_value_array[id_sequence] = p_value
            if biosensor == 1:
                biosensor_array[id_sequence] = 1
                print(id_sequence)

non_biosensor_array = biosensor_array == 0
biosensor_array = biosensor_array == 1

# fold_change_array = np.array(
#     [sequence_info[1] for sequence_info in uniq_sequences_info], dtype=np.float64)
# fold_change_sd_array = np.array(
#     [sequence_info[2] for sequence_info in uniq_sequences_info], dtype=np.float64)
# p_value_array = np.array(
#     [sequence_info[3] for sequence_info in uniq_sequences_info], dtype=np.float64)


p_value_minus_log10_array = - np.log10(p_value_array)

# fig1, ax = plt.subplots()
# # ax.scatter(fold_change_array, fold_change_sd_array, s=2, alpha=1)
# ax.scatter(fc_array[non_biosensor_array],
#            p_value_minus_log_array[non_biosensor_array], s=2, alpha=1)
# ax.scatter(fc_array[biosensor_array],
#            p_value_minus_log_array[biosensor_array], s=2, alpha=1, color="red")
# ax.grid()

# ax.set_xlabel("fold change")
# ax.set_ylabel("-log10(p-value)")
# # ax.set_xlim([0.2, 10])
# ax.set_ylim([0.0, 2])
# ax.set_xscale('log')
# ax.set_yscale('log')
# fig.savefig("fold_change_vs_fold_change_se.pdf", dpi=300)


fig2, ax2 = plt.subplots()
ax2.scatter(fc_array[non_biosensor_array],
            fc_se_array[non_biosensor_array], s=2, alpha=1, color="#057D54")
ax2.scatter(fc_array[biosensor_array],
            fc_se_array[biosensor_array], s=2, alpha=1, color="#9B0138")

ax2.set_xlabel("Fold change")
ax2.set_ylabel("Standard error")


fig2.savefig(f"{snakemake.output[0]}", format="svg", dpi=1200)


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

fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
ax1.scatter(fraction_negative[non_biosensor_array],
            fraction_positive[non_biosensor_array], s=5, alpha=1, color="#057D54")
ax1.scatter(fraction_negative[biosensor_array],
            fraction_positive[biosensor_array], s=5, alpha=1, color="red")

ax1.plot(np.arange(100), np.arange(100),
         color="#9B0138", label="reference", linestyle="dashed", linewidth=1)


ax1.plot(np.arange(100), np.arange(100),
         color="#9B0138", label="reference", linestyle="dashed", linewidth=1)
ax1.set_xlim([0, 100])
ax1.set_ylim([0, 100])
ax1.xaxis.set_major_locator(MultipleLocator(10))
ax1.xaxis.set_minor_locator(MultipleLocator(5))
ax1.yaxis.set_major_locator(MultipleLocator(10))
ax1.yaxis.set_minor_locator(MultipleLocator(5))

# Add a percentage sign to the values
ax1.xaxis.set_major_formatter(FormatStrFormatter("%d%%"))
ax1.yaxis.set_major_formatter(FormatStrFormatter("%d%%"))

ax1.set_xlabel("Fraction cleaved (-ligand)")
ax1.set_ylabel("Fraction cleaved (+ligand)")


fig1.savefig(f"{snakemake.output[1]}", format="svg", dpi=1200)
