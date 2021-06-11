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
    cs_sd = calc_helpers.calculate_sample_standard_deviation(cs_pos)

    # Store the 5% and 95% percentiles of the mean cleavage fraction
    cs_5, cs_95 = calc_helpers(np.mean(cs), cs_sd)
    return cs_sd, cs_5, cs_95


def bootstrap_fold_change_with_replacement(samples_pos, samples_neg, num_samples=1000):
    """
    Function does bootstrapping on the data.
    args:\n
    samples_pos: (array) cleaved prefix array yes/no of reads with condition ligand present.\n
    \texample: let's say we have a sequence where three prefixes were cleaved and two uncleaved
    then the array will look as follows: [0, 0, 1, 1, 1].\n
    samples_neg: (array) cleaved prefix array yes/no of reads with condition ligand not present.\n
    num_samples: (int) number of bootstrap samples to take.\n
    \n
    returns:

    """
    n = samples.shape[0]  # gives size of the array

    # First do the positive samples
    # Create a mask which indexes to bootstrap with replacement
    random_mask_pos = random.randint(0, high=n, size=(
        num_samples, n), dtype=np.int32)    # high value is exclusive

    bootstrap_samples_pos = samples[random_mask_pos]
    # Get number of prefixes
    num_clvd_s = np.count_nonzero(bootstrap_samples_pos == 1, axis=1)
    num_unclvd_s = np.count_nonzero(bootstrap_samples_pos == 0, axis=1)

    cs_pos = calc_helpers.calculate_cleavage_fraction(
        num_clvd_s, r_clvd_ref, num_unclvd_s, r_unclvd_ref)

    # Calculate standard deviation
    cs_pos_sd = calc_helpers.calculate_sample_standard_deviation(cs_pos)

    # Store the 5% and 95% percentiles of the mean cleavage fraction
    cs_pos_5, cs_pos_95 = calc_helpers(np.mean(cs_pos), cs_pos_sd)

    # Then do the negative samples
    # Create a mask which indexes to bootstrap with replacement
    random_mask_neg = random.randint(0, high=n, size=(
        num_samples, n), dtype=np.int32)    # high value is exclusive

    bootstrap_samples_neg = samples[random_mask_neg]
    # Get number of prefixes
    num_clvd_s = np.count_nonzero(bootstrap_samples_neg == 1, axis=1)
    num_unclvd_s = np.count_nonzero(bootstrap_samples_neg == 0, axis=1)

    cs_neg = calc_helpers.calculate_cleavage_fraction(
        num_clvd_s, r_clvd_ref, num_unclvd_s, r_unclvd_ref)

    bootstrap_means = np.mean(bootstrap_samples, axis=1)
