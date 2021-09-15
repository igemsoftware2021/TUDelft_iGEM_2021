import numpy as np


def calc_sample_standard_deviation(samples):
    """Function calculates the standard deviation of a sample. Formula for a sample standard
    deviation:\n
    \n
    s = sqrt(sum(x-x_mean)**2)/(n-1)).\n
    """
    # s = sqrt(sum((x-x_mean)**2)/n-1)
    n = samples.shape[0]
    sample_mean = np.mean(samples)
    sample_standard_deviation = np.sqrt(
        np.sum(((samples - sample_mean)**2))/(n-1))
    return sample_standard_deviation


def calc_95_confidence_interval(mean, standard_deviation):
    """Function returns a 95% confidence interval."""
    return (mean - 1.96 * standard_deviation, mean + 1.96 * standard_deviation)


def calc_cleavage_fraction(r_clvd_s: int, r_clvd_ref: int, r_unclvd_s: int, r_unclvd_ref: int):
    """Function calculates the cleavage fraction. The formula is from the DRIVER paper page 12, equation (1).\n
    args:\n
    r_clvd_s: (int) number of reads of sequence with cleaved prefix.\n
    r_clvd_ref: (int) number of reads of reference sequences with cleaved prefix.\n
    r_unclvd_s: (int) number of reads of sequence with uncleaved prefix.\n
    r_unclvd_ref: (int) number of reads of reference sequences with uncleaved prefix.\n
    \n
    returns:\n
    (float) cleavage fraction value
    """
    return (r_clvd_s/r_clvd_ref) / ((r_unclvd_s/r_unclvd_ref) + (r_clvd_s/r_clvd_ref))


def calc_fold_change(cs_pos: float, cs_neg: float, k: float = 1.0):
    """Function calculates the fold_change. The formula is from the DRIVER paper page 12, equation (2).\n
    args:\n
    cs_pos: (float) cleavage fraction of sequence under condition ligand present.\n
    cs_neg: (float) cleavage fraction of sequence under condition ligand not present.\n
    k: (int, float) factor used to compansate for slight variations in the experimental conditions.\n
    \n
    returns:\n
    (float) fold change value
    """
    return k * ((1 - cs_pos)/(1 - cs_neg))