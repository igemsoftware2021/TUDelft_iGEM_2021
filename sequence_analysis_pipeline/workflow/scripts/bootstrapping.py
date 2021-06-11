import numpy as np
import calc_helpers


def bootstrap_cleavage_fraction_with_replacement(data, r_clvd_ref, r_unclvd_ref, num_samples=1000):

    # Create a mask which indexes to bootstrap with replacement
    random_mask = random.randint(0, high=n, size=(
        num_samples, n), dtype=np.int32)    # high value is exclusive

    bootstrap_samples = data[random_mask]
    # Get number of prefixes
    r_clvd_s = np.count_nonzero(bootstrap_samples == 1, axis=1)
    r_unclvd_s = np.count_nonzero(bootstrap_samples == 0, axis=1)

    cs = calc_helpers.calculate_cleavage_fraction(
        r_clvd_s, r_clvd_ref, r_unclvd_s, r_unclvd_ref)

    # Calculate standard deviation of the bootstrapped cleavage fractions
    cs_sd = calc_helpers.calculate_sample_standard_deviation(cs)

    # Store the 5% and 95% percentiles of the mean cleavage fraction
    cs_5_perc, cs_95_perc = calc_helpers(np.mean(cs), cs_sd)
    return cs, cs_sd, cs_5_perc, cs_95_perc


def bootstrap_fold_change_with_replacement(cs_pos, cs_neg, r_clvd_ref, r_unclvd_ref, k=1, num_samples=1000):
    """
    Function does bootstrapping for cleavage fraction and fold change on the data.
    args:\n
    samples_pos: (array) cleaved prefix array yes/no of reads with condition ligand present.\n
    \texample: let's say we have a sequence where three prefixes were cleaved and two uncleaved
    then the array will look as follows: [0, 0, 1, 1, 1].\n
    samples_neg: (array) cleaved prefix array yes/no of reads with condition ligand not present.\n
    num_samples: (int) number of bootstrap samples to take.\n
    \n
    returns:

    """
    if cs_pos.shape != cs_neg.shape:
        raise RuntimeError(
            "The shapes of cs_pos array and cs_neg array are not the same.")

    # TODO decide on having the bootstrap_cleavage_fraction_with_replacement() in the function yes or no.
    # Positive data first
    cs_pos, cs_pos_sd, cs_pos_5_perc, cs_pos_95_perc = bootstrap_cleavage_fraction_with_replacement(
        data_pos, r_clvd_ref, r_unclvd_ref, num_samples=num_samples)

    # Then negative data
    cs_neg, cs_neg_sd, cs_neg_5_perc, cs_neg_95_perc = bootstrap_cleavage_fraction_with_replacement(
        data_neg, r_clvd_ref, r_unclvd_ref, num_samples=num_samples)

    fold_changes = calc_helpers.calculate_fold_change(cs_pos, cs_neg, k=k)
    fold_change_sd = calc_helpers.calculate_sample_standard_deviation(
        fold_changes)
    fold_change_5_perc, fold_change_95_perc = calc_helpers.calculate_95_confidence_interval(
        np.mean(fold_changes), fold_change_sd)

    return fold_change_sd, fold_change_5_perc, fold_change_95_perc
