import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Folate
# ------
# ------
#
# Reference sequences
# -------------------
# +ligand, uncleaved = 586
# +ligand, cleaved = 494
# -ligand, uncleaved = 399
# -ligand, cleaved = 333

# Sequence 1
# ----------
# +ligand, uncleaved = 301
# +ligand, cleaved = 52
# -ligand, uncleaved = 185
# -ligand, cleaved = 1349
#
# Sequence 5
# ----------
# +ligand, uncleaved = 2007
# +ligand, cleaved = 35
# -ligand, uncleaved = 377
# -ligand, cleaved = 1182
#
#
#
# Vitamins
# ------
# ------
#
# Reference sequences
# -------------------
# +ligand, uncleaved = 819
# +ligand, cleaved = 624
# -ligand, uncleaved = 627
# -ligand, cleaved = 457

# Sequence 0
# ----------
# +ligand, uncleaved = 603
# +ligand, cleaved = 67
# -ligand, uncleaved = 382
# -ligand, cleaved = 1786
#
# Sequence 2
# ----------
# +ligand, uncleaved = 13118
# +ligand, cleaved = 330
# -ligand, uncleaved = 2678
# -ligand, cleaved = 9622
#
# Sequence 3
# ----------
# +ligand, uncleaved = 2298
# +ligand, cleaved = 44
# -ligand, uncleaved = 295
# -ligand, cleaved = 1021


folate_r_ref_unclvd_pos = 586
folate_r_ref_clvd_pos = 494
folate_r_ref_unclvd_neg = 399
folate_r_ref_clvd_neg = 333

vitamins_r_ref_unclvd_pos = 819
vitamins_r_ref_clvd_pos = 624
vitamins_r_ref_unclvd_neg = 627
vitamins_r_ref_clvd_neg = 457

# Structure for one sequence
# [(+ligand, uncleaved), (+ligand, cleaved), (-ligand, uncleaved), (-ligand, cleaved)]
r_seqs = [[301, 52, 185, 1349], [2007, 35, 377, 1182], [
    603, 67, 382, 1786], [13118, 330, 2678, 9622], [2298, 44, 295, 1021]]

seq_names = ["Sequence A", "Sequence B",
             "Sequence C", "Sequence D", "Sequence E"]

fig = plt.figure(figsize=(12, 8), dpi=125)
ax1 = plt.subplot2grid(shape=(2, 14), loc=(0, 0), colspan=4)
ax6 = plt.subplot2grid(shape=(2, 14), loc=(0, 4), colspan=1)
ax2 = plt.subplot2grid(shape=(2, 14), loc=(0, 5), colspan=4)
ax7 = plt.subplot2grid(shape=(2, 14), loc=(0, 9), colspan=1)
ax3 = plt.subplot2grid(shape=(2, 14), loc=(0, 10), colspan=4)
ax4 = plt.subplot2grid(shape=(2, 14), loc=(1, 2), colspan=4)
ax8 = plt.subplot2grid(shape=(2, 14), loc=(1, 6), colspan=2)
ax5 = plt.subplot2grid(shape=(2, 14), loc=(1, 8), colspan=4)


ax6.remove()
ax7.remove()
ax8.remove()

handles = [ax1, ax2, ax3, ax4, ax5]


# Width of a bar in a bar plot
width = 0.025

max_ind = 0.1
ind = np.array([0, max_ind])

for i in range(5):

    # Plot handle
    handle = handles[i]

    # Sequence read counts
    r_seq = r_seqs[i]
    seq_name = seq_names[i]

    if i < 2:
        r_ref_unclvd_pos = folate_r_ref_unclvd_pos
        r_ref_clvd_pos = folate_r_ref_clvd_pos
        r_ref_unclvd_neg = folate_r_ref_unclvd_neg
        r_ref_clvd_neg = folate_r_ref_clvd_neg
    else:
        r_ref_unclvd_pos = vitamins_r_ref_unclvd_pos
        r_ref_clvd_pos = vitamins_r_ref_clvd_pos
        r_ref_unclvd_neg = vitamins_r_ref_unclvd_neg
        r_ref_clvd_neg = vitamins_r_ref_clvd_neg

    # Plot uncleaved
    handle.bar(ind-0.55*width, [r_seq[0]/r_ref_unclvd_pos, r_seq[2] /
                                r_ref_unclvd_neg], width=width, color="#9B0138", label="Uncleaved")
    # Plot cleaved
    handle.bar(ind+0.55*width, [r_seq[1]/r_ref_clvd_pos, r_seq[3] /
                                r_ref_clvd_neg], width=width, color="#057D54", label="Cleaved")
    handle.set_xticks([0, max_ind])
    handle.set_xticklabels(["+ligand", "-ligand"])
    handle_ylim = handle.get_ylim()
    handle.set_ylim((0, handle_ylim[1]*1.25))

    if i != 3:
        handle.yaxis.set_major_locator(MultipleLocator(1))
        handle.yaxis.set_minor_locator(MultipleLocator(0.1))
    else:
        handle.yaxis.set_major_locator(MultipleLocator(5))
        handle.yaxis.set_minor_locator(MultipleLocator(1))

    handle.set_ylabel("$r_{\mathrm{s}}\;/\;r_{\mathrm{ref}}$")
    handle.set_xlabel(f"{seq_name}", fontweight="bold")
    handle.legend(loc="upper left")

# colors = {"uncleaved": "#9B0138", "cleaved": "#057D54"}
# labels = list(colors.keys())
# legend_handles = [plt.Rectangle((0, 0), 1, 1, color=colors[label])
#                   for label in labels]

# fig.legend(legend_handles, labels, loc="upper center")

plt.show()
