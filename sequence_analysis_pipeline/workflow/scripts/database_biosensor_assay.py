import scipy.stats
from database_interface import DatabaseInterfaceCleanSequences

# We need to do one-sided

seq_fold_change = 3.0
expected_fold_change = 1
z_score = (seq_fold_change - expected_fold_change) / 1.5
print(z_score)
p_value = scipy.stats.norm.sf(abs(z_score))  # one-sided
print(p_value)
