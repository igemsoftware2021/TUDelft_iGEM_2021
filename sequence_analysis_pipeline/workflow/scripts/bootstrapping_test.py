import numpy as np

samples = np.zeros(100)
samples[-25:] = 1
# samples = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
highest_index = samples.shape[0]   # len gives highest_index + 1
print(highest_index)
mask = np.random.randint(0, high=highest_index, size=(
    10, highest_index), dtype=np.int32)    # high value is exclusive

new_values = samples[mask]
num_clvd_s = np.count_nonzero(new_values == 1, axis=1)
print(num_clvd_s)
# print(mask)
# print(new_values)

bootstrap_means = np.mean(new_values, axis=1)
# print(bootstrap_means)
n = bootstrap_means.shape[0]
sample_mean = np.mean(samples)
print(bootstrap_means, sample_mean)
print("test", np.sum((bootstrap_means - sample_mean)**2))
ssd = np.sqrt(np.sum((bootstrap_means - sample_mean)**2)/(n-1))
print(ssd)
sdd_sqrt_n = ssd / np.sqrt(n)
print(sdd_sqrt_n)
confidence_interval = (sample_mean - 1.96*ssd, sample_mean + 1.96*ssd)
print(confidence_interval)


test_a = np.zeros(20, dtype=np.int8)
test_a[-10:] = 1
print(test_a)
