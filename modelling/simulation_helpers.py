import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from models import model_prokaryotic_readout
from alignment_helpers import weighted_running_average
from standard_values import standard_parameters_prokaryotic, standard_constants


def create_absorbance_sd(num_points, low=0.0, high=2.0, linear_coef=0.01, random_coef=0.0):
    """Function creates an array with absorbance values and an array with the corresponding standard deviation"""
    absorbance = np.linspace(low, high, num=num_points, dtype=np.float32)
    sd = linear_coef * absorbance + random_coef * np.random.randn(num_points)
    return absorbance, sd


def find_nearest_idx(array, value):
    """Function finds index of the array where the value is closest to the given value.

    Example:
    array = np.array([0.0, 0.6, 0.9, 1.0, 1.5, 5.0, 7.6])
    value = 1.4

    returns idx = 4
    """
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array-value))
    return idx


def find_nearest(array, value):
    """Returns the value in the array that is closest to the given value

    Example:
    array = np.array([0.0, 0.6, 0.9, 1.0, 1.5, 5.0, 7.6])
    value = 1.4

    returns 1.5, since value at index = 4 is the closest
    """
    idx = find_nearest_idx(array, value)
    return array[idx]


def simulate_time_lag(time, absorbance_data, t_lag=10.0, t_sd_lag=5.0):
    """Function simulates time lag in the cell free systems"""
    time_lag = t_lag + t_sd_lag * np.random.randn()
    dt = time[1] - time[0]

    time_lag_array = np.arange(0.0, time_lag, step=dt, dtype=np.float32)
    absorbance_lag_array = np.zeros(time_lag_array.shape, dtype=np.float32)

    time = np.concatenate((time_lag_array, time), axis=0)
    absorbance_data = np.concatenate(
        (absorbance_lag_array, absorbance_data), axis=0)
    return time, absorbance_data


def simulate_time_measurement(time, absorbance_data, dt=1.0, dt_sd=0.1):
    # First fit a spline to the data
    spl = UnivariateSpline(time, absorbance_data, k=1, s=0, ext=2)

    # Create the new time array
    # First find the total measuring time
    time_total = time[-1]

    # Create a list where we will store the new timepoints
    time_measurement = []

    # Create a temporary variable to store how far in time you are
    time_temp = 0
    while time_temp < time_total:
        time_measurement.append(time_temp)
        time_temp = time_temp + (dt + dt_sd * np.random.randn())

    # Transform the list to a numpy array
    time_measurement = np.array(time_measurement, dtype=np.float32)
    absorbance_data_time_measurement = spl(time_measurement)

    return time_measurement, absorbance_data_time_measurement


def simulate_absorbance_noise(absorbance_data, absorbance_value, absorbance_sd):
    """Function generates absorbance noise absorbance_data is the actual array with the real data
    absorbance_value is an array with absorbance values and absorbance_sd gives a standard deviation for
    every absorbance value"""
    absorbance_data_noise = np.copy(absorbance_data)

    # Add noise depending on the standard deviations for different absorbance values
    for i in range(0, absorbance_data_noise.shape[0]):
        idx = find_nearest_idx(absorbance_value, absorbance_data_noise[i])
        absorbance_data_noise[i] = absorbance_data_noise[i] + \
            absorbance_sd[idx] * np.random.randn()

    return absorbance_data_noise


def simulate_hardware(time, absorbance_data, absorbance_value, absorbance_sd, dt: float = 1.0, dt_sd: float = 0.1, t_lag: float = None, t_sd_lag: float = None):
    """Function simulates the hardware, so it removes time points and it adds noise to the absorbance
    Only time lag is added if t_lag and t_sd_lag are not None
    """

    if t_lag is not None and t_sd_lag is not None:
        time, absorbance_data = simulate_time_lag(
            time, absorbance_data, t_lag=t_lag, t_sd_lag=t_sd_lag)

    time_noise, absorbance_data = simulate_time_measurement(
        time, absorbance_data, dt=dt, dt_sd=dt_sd)
    absorbance_data_noise = simulate_absorbance_noise(
        absorbance_data, absorbance_value, absorbance_sd)

    return time_noise, absorbance_data_noise


def simulate_noise(time, data, absorbance_values, sd_values, dt=1, dt_sd=0.1):
    # First fit a spline to the data
    f = UnivariateSpline(time, data, k=1, s=0, ext=2)

    time_new2 = np.arange(0, time[-1], step=0.1)
    data_new2 = f(time_new2)

    # Create the figure and the line
    fig1, ax1 = plt.subplots()
    ax1.plot(time, data)
    ax1.plot(time_new2, data_new2)
    fig1.show()
    plt.show()

    # First determine the total time
    total_time = time[-1]

    temp_time = 0
    time_new = [0]
    while temp_time < total_time:
        time_to_add = dt + dt_sd * np.random.randn()
        temp_time = temp_time + time_to_add
        if temp_time < total_time:
            time_new.append(temp_time)

    time_new = np.array(time_new, dtype=np.float32)
    data_new = f(time_new)

    # Add noise depending on the standard deviations for different absorbance values
    for i in range(0, data_new.shape[0]):
        idx = find_nearest_idx(absorbance_values, data_new[i])
        data_new[i] = data_new[i] + sd_values[idx] * np.random.randn()

    return time_new, data_new


if __name__ == "__main__":
    # Parameters
    parameters = standard_parameters_prokaryotic()
    # Constants
    constants = standard_constants()

    t_tot = 7200  # total time [s]
    dt = 0.01  # timestep [s]

    # Initial concentration of the beta-galactosidase gene [μM]
    dna_i = 5*10**-3
    s_i = 150         # Initial substrate concentration [μM]

    initial_conditions = np.array([dna_i, s_i, 500*10**-3])

    # Running the model
    time, absorbance = model_prokaryotic_readout(
        parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

    absorbance_values, sd_values = create_absorbance_sd(30)

    f = UnivariateSpline(absorbance_values, sd_values, k=3)

    x_new = np.arange(0, 2.0, step=0.01, dtype=np.float32)
    y_new = f(x_new)

    time_new, absorbance_new = simulate_hardware(
        time, absorbance, absorbance_values, sd_values)

    time_weighted, absorbance_weighted = weighted_running_average(
        time_new, absorbance_new, n=121)
    # Create the figure and the line
    fig1, ax1 = plt.subplots()
    ax1.plot(time, absorbance)
    ax1.scatter(time_new, absorbance_new, s=1, alpha=0.5, color="red")
    ax1.plot(time_weighted, absorbance_weighted, color="green")
    fig1.show()
    plt.show()
