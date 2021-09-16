import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from models import model_prokaryotic_readout
from alignment_helpers import weighted_running_average


def create_absorbance_sd(num_points, low=0.0, high=2.0, linear_coef=0.005, random_coef=0.002):
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
    # (#) denotes the position in the parameters array
    k_ts = 6.7*10**-5     # (0)  Transcription rate [uM/s]
    k_tl = 7.2*10**-5     # (1)  Enzyme translation rate [uM/s]
    k_mat = 0.5*10**-3    # (2)  Maturation rate of beta-galactosidase [1/s]
    k_cat = 5.14*10**1    # (3)  Catalytic rate of beta-galactosidase [1/s]
    k_s = 8.5*10**-3      # (4)  Michaelis constant of transcription [μM]
    # (5)  Scaling factor for the transcription resources [1/s]
    kc_s = 1.8*10**-4
    k_l = 65.8*10**-3     # (6)  Michaelis constant of translation [μM]
    k_tlr = 6*10**-5      # (7) Michaelis constant of translation resources [-]
    k_m = 50              # (8) Michaelis constant of beta-galactosidase [μM]
    deg_mrna = 1.3*10**-5  # (9) Degradation rate of mRNA [1/s]
    # (10) Degradation rate of translation resources [1/s]
    deg_tlr = 7.5*10**-5
    k_on = 1*10**-2       # (11)  Association rate of vitamin and umRNA [1/µMs]
    k_off = 1*10**-2      # (12)  Dissociation rate of vitamin and umRNA [1/s]
    k_c = (1/60)/10            # (13)  Cleaving rate of umRNA [1/s]
    parameters = np.array([k_ts, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                           k_tlr, k_m, deg_mrna, deg_tlr, k_on, k_off, k_c])  # Array containing above parameters

    t_tot = 7200  # total time [s]
    dt = 0.01  # timestep [s]

    # Constants
    # (#) denotes the position in the constants array
    h = 0.020  # (0) Height of the paper [cm]
    # (1) Exctinction coefficient of CPRG at a wavelength of 410 in 25% serum [1/(µM*cm)]
    eps_cprg = 0.294
    # (2) Exctinction coefficient of CPR at a wavelength of 580 in 25% serum [1/(µM*cm)]
    eps_cpr = 0.539
    # I_null = ... Has to be filled in after the experiments.

    # Array containing above constants
    constants = np.array([h, eps_cprg, eps_cpr])

    # Initial concentration of the beta-galactosidase gene [μM]
    dna_i = 5*10**-3
    s_i = 150         # Initial substrate concentration [μM]

    initial_conditions = np.array([dna_i, s_i, 500*10**-3])

    # Running the model
    time, absorbance = model_prokaryotic_readout(
        parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

    absorbance_values, sd_values = absorbance_sd(30)

    f = UnivariateSpline(absorbance_values, sd_values, k=3)

    x_new = np.arange(0, 2.0, step=0.01, dtype=np.float32)
    y_new = f(x_new)

    time_new, absorbance_new = simulate_noise(
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
