from helpers import read_from_csv
import matplotlib.pyplot as plt

path = "hardware/test_results/"
timepoints_plot, absorbance_plot = read_from_csv(path, "sensor_1")
fig1, ax1 = plt.subplots()
ax1.plot(timepoints_plot, absorbance_plot)
fig1.show()

timepoints_plot, temperature_error = read_from_csv(path, "temperature_error")
fig2, ax2 = plt.subplots()
ax2.plot(timepoints_plot, temperature_error)
fig2.show()
print(timepoints_plot)

plt.show()
