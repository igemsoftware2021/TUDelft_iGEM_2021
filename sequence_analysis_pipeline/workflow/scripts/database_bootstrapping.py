

# Positive data first
cs_pos, cs_pos_sd, cs_pos_5_perc, cs_pos_95_perc = bootstrap_cleavage_fraction_with_replacement(
    data_pos, r_clvd_ref, r_unclvd_ref, num_samples=num_samples)

# Then negative data
cs_neg, cs_neg_sd, cs_neg_5_perc, cs_neg_95_perc = bootstrap_cleavage_fraction_with_replacement(
    data_neg, r_clvd_ref, r_unclvd_ref, num_samples=num_samples)
