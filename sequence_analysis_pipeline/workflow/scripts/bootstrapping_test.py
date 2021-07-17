from random import sample
import numpy as np
from calc_helpers import calc_cleavage_fraction
from bootstrapping import bootstrap_cleavage_fraction_with_replacement, bootstrap_fold_change_with_replacement
import matplotlib.pyplot as plt

r_unclvd_s_neg = 320
r_clvd_s_neg = 3201

r_clvd_ref = 100
r_unclvd_ref = 120


r_unclvd_s_pos = 3501
r_clvd_s_pos = 432

cs_neg_mean, cs_neg_sd, cs_pos_mean, cs_pos_sd, fold_changes, fold_change_mean, fold_change_sd = bootstrap_fold_change_with_replacement(
    r_clvd_s_neg, r_unclvd_s_neg, r_clvd_s_pos, r_unclvd_s_pos, r_clvd_ref, r_unclvd_ref, k=0.55, num_samples=1000)

print(cs_neg_mean, cs_neg_sd, cs_pos_mean, cs_pos_sd)
print(fold_changes)
print(fold_change_mean)
print(fold_change_sd)

plt.hist(fold_changes, bins=20)
plt.show()
