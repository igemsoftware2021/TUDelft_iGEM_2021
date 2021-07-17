import numpy as np
import calc_helpers


# def bootstrap_cleavage_fraction_with_replacement(data, r_clvd_ref, r_unclvd_ref, num_samples=1000):

#     n = data.shape[0]
#     # Create a mask which indexes to bootstrap with replacement
#     random_mask = np.random.randint(0, high=n, size=(
#         num_samples, n), dtype=np.int32)    # high value is exclusive

#     bootstrap_samples = data[random_mask]
#     # Get number of prefixes
#     r_clvd_s = np.count_nonzero(bootstrap_samples == 1, axis=1)
#     r_unclvd_s = np.count_nonzero(bootstrap_samples == 0, axis=1)

#     print(r_clvd_s)

#     cs = calc_helpers.calc_cleavage_fraction(
#         r_clvd_s, r_clvd_ref, r_unclvd_s, r_unclvd_ref)

#     # Calculate estimated mean
#     cs_mean = np.mean(cs)

#     # Calculate standard deviation of the bootstrapped cleavage fractions
#     cs_sd = calc_helpers.calc_sample_standard_deviation(cs)

#     return cs, cs_mean, cs_sd


def bootstrap_fold_change_with_replacement(r_clvd_s_neg, r_unclvd_s_neg, r_clvd_s_pos, r_unclvd_s_pos, r_clvd_ref_neg, r_unclvd_ref_neg, r_clvd_ref_pos, r_unclvd_ref_pos, k=1, num_samples=1000) -> tuple:
    # TODO update doc string
    """
    Function does bootstrapping for cleavage fraction and fold change on the data.
    args:\n
    samples_pos: (array) cleaved prefix array yes/no of reads with condition ligand present.\n
    \texample: let's say we have a sequence where three prefixes were cleaved and two uncleaved
    then the array will look as follows: [0, 0, 1, 1, 1].\n
    samples_neg: (array) cleaved prefix array yes/no of reads with condition ligand not present.\n
    num_samples: (int) number of bootstrap samples to take.\n
    \n
    returns:\n
    A tuple containing the following values in the same order:\n
    cs_pos_sd: (float) sample standard deviaton of the cleavage fraction for the condition ligand present.\n
    cs_pos_5_perc: (float) lower 95% limit for the cleavage fraction for the condition ligand present.\n
    cs_pos_95_perc: (float) upper 95% limit for the cleavage fraction for the condition ligand present.\n
    cs_neg_sd: (float) sample standard deviation of the cleavage fraction for the condition ligand not present.\n
    cs_neg_5_perc: (float) lower 95% limit for the cleavage fraction for the condition ligand not present.\n
    cs_neg_95_perc: (float) upper 95% limit for the cleavage fraction for the condition ligand not present.\n
    fold_change_sd: (float) sample standard deviation of the fold change.\n
    fold_change_se: (float) standard error of the fold change.\n
    fold_change_5_perc: (float) lower 95% limit for the fold change.\n
    fold_change_95_perc: (float) upper 95% limit for the fold change.\n
    """

    # Create an array with all the read counts
    # A number is linked to the following read:
    # 0: r_clvd_s_neg
    # 1: r_unclvd_s_neg
    # 2: r_clvd_s_pos
    # 3: r_unclvd_s_pos
    sample_data = np.zeros(
        r_clvd_s_neg + r_unclvd_s_neg + r_clvd_s_pos + r_unclvd_s_pos, dtype=np.int8)
    # Fill the array with an amount of r_seq_clvd_neg ones.
    sample_data[r_clvd_s_neg:(r_clvd_s_neg + r_unclvd_s_neg)] = 1
    sample_data[(r_clvd_s_neg + r_unclvd_s_neg):(r_clvd_s_neg + r_unclvd_s_neg + r_clvd_s_pos)] = 2
    sample_data[-r_unclvd_s_pos:] = 3

    n = sample_data.shape[0]
    # Create a mask that has random idx to bootstrap with replacement
    random_idx = np.random.randint(0, high=n, size=(
        num_samples, n), dtype=np.int32)    # high value is exclusive

    bootstrap_samples = sample_data[random_idx]

    # Get read_counts for every prefix and condition
    sample_r_clvd_s_neg = np.count_nonzero(bootstrap_samples == 0, axis=1)
    sample_r_unclvd_s_neg = np.count_nonzero(bootstrap_samples == 1, axis=1)
    sample_r_clvd_s_pos = np.count_nonzero(bootstrap_samples == 2, axis=1)
    sample_r_unclvd_s_pos = np.count_nonzero(bootstrap_samples == 3, axis=1)

    cs_neg = calc_helpers.calc_cleavage_fraction(
        sample_r_clvd_s_neg, r_clvd_ref_neg, sample_r_unclvd_s_neg, r_unclvd_ref_neg)

    # Calculate estimated mean
    cs_neg_mean = np.mean(cs_neg)

    # Calculate standard deviation of the bootstrapped cleavage fractions
    cs_neg_se = calc_helpers.calc_sample_standard_deviation(cs_neg)

    cs_pos = calc_helpers.calc_cleavage_fraction(
        sample_r_clvd_s_pos, r_clvd_ref_pos, sample_r_unclvd_s_pos, r_unclvd_ref_pos)

    # Calculate estimated mean
    cs_pos_mean = np.mean(cs_pos)

    # Calculate standard deviation of the bootstrapped cleavage fractions
    cs_pos_se = calc_helpers.calc_sample_standard_deviation(cs_pos)

    fold_changes = calc_helpers.calc_fold_change(cs_pos, cs_neg, k=k)

    fold_change_mean = np.mean(fold_changes)

    fold_change_se = calc_helpers.calc_sample_standard_deviation(
        fold_changes)

    return cs_neg_mean, cs_neg_se, cs_pos_mean, cs_pos_se, fold_changes, fold_change_mean, fold_change_se
