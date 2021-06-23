import numpy as np
import calc_helpers


def bootstrap_cleavage_fraction_with_replacement(data, r_clvd_ref, r_unclvd_ref, num_samples=1000):

    n = data.shape[0]
    # Create a mask which indexes to bootstrap with replacement
    random_mask = np.random.randint(0, high=n, size=(
        num_samples, n), dtype=np.int32)    # high value is exclusive

    bootstrap_samples = data[random_mask]
    # Get number of prefixes
    r_clvd_s = np.count_nonzero(bootstrap_samples == 1, axis=1)
    r_unclvd_s = np.count_nonzero(bootstrap_samples == 0, axis=1)

    cs = calc_helpers.calc_cleavage_fraction(
        r_clvd_s, r_clvd_ref, r_unclvd_s, r_unclvd_ref)

    # Calculate standard deviation of the bootstrapped cleavage fractions
    cs_sd = calc_helpers.calc_sample_standard_deviation(cs)

    # Store the 5% and 95% percentiles of the mean cleavage fraction
    cs_5_perc, cs_95_perc = calc_helpers(np.mean(cs), cs_sd)
    return cs, cs_sd, cs_5_perc, cs_95_perc


def bootstrap_fold_change_with_replacement(cs_neg, cs_pos, k=1) -> tuple:
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

    if cs_neg.shape != cs_pos.shape:
        raise ValueError(
            "Shape of cs_neg is not the same as shape cs_pos, check if the function bootstrap_fold_change_with_replacement has the correct input values.")

    num_samples = cs_neg.shape[0]

    fold_changes = calc_helpers.calc_fold_change(cs_pos, cs_neg, k=k)
    fold_change_sd = calc_helpers.calc_sample_standard_deviation(
        fold_changes)

    fold_change_5_perc, fold_change_95_perc = calc_helpers.calc_95_confidence_interval(
        np.mean(fold_changes), fold_change_sd)

    # fold change standard error
    fold_change_se = fold_change_sd / np.sqrt(num_samples)

    return fold_change_sd, fold_change_se, fold_change_5_perc, fold_change_95_perc
