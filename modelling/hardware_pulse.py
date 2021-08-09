from re import X
import numpy as np
import scipy.signal
import scipy.special
import scipy.fft
import matplotlib.pyplot as plt

cutoff_freq = 3.3 * 10**6

# Frequency
freq = 500 * 10**3

# Sample spacing
T = 1.0/freq

# Sample rate
SAMPLE_RATE = 100 * freq
# Duration
DURATION = 0.01  # s

# Number of sample points
N = int(DURATION * SAMPLE_RATE)
print(N)
x = np.linspace(0, DURATION, N, endpoint=False, dtype=np.float32)
y = scipy.signal.square(2 * np.pi * freq * x)

y = np.maximum(y, np.zeros(y.shape, dtype=np.float32))

yf = scipy.fft.rfft(y)
xf = scipy.fft.fftfreq(N, 1/SAMPLE_RATE)[:N//2]

xf_temp = np.abs(xf - cutoff_freq)
idx = np.argmin(xf_temp)
print(idx)
yf_adapted = np.copy(yf)
yf_adapted[idx:] = 0.0
yf_adapted[((N//2)+idx):] = 0.0

y_adapted = scipy.fft.irfft(yf_adapted)

# plt.plot(xf, 2.0/N * np.abs(yf_adapted))
plt.plot(x[:int(T*SAMPLE_RATE*5)], y_adapted[:int(T*SAMPLE_RATE*5)])
plt.show()
