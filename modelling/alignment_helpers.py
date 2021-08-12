import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d


def running_average(a, n=11, mode="same"):
    return np.convolve(a, np.ones(n)/n, mode=mode)


def finite_difference(x, y):
    """
    Compute the finite difference between y values.

    Parameters
    ----------
    x: (array) this array is used to determine the delta_x
    y: (array) finite difference is calculated for all points

    Returns
    -------
    (np.array) all the computed finite differences.
    """

    # Calculate the central difference dydx
    delta_y = np.convolve(y, np.array(
        [1, 0, -1], dtype=np.float64), mode="valid")
    # Calculate the delta x for each datapoint, because the time between measurements can change.
    delta_x = np.convolve(x, np.array(
        [1, 0, -1], dtype=np.float64), mode="valid")
    dydx = np.divide(delta_y, delta_x)

    return dydx


def find_platau_index(x, y, tol=1e-8, n=11):

    # Calculate the dydx for the graph
    dydx = finite_difference(x, y)

    # Calculate the running average for the found dydx
    dydx_average = running_average(dydx, n=n)

    start_idx = dydx_average.shape[0]//2

    for i in range(start_idx, dydx_average.shape[0]):
        if dydx_average[i] < tol:
            return i + 2  # this is not the right one because the array got smaller, +2 ?


def find_start_index(x, y, tol=1e-8, n=11):

    # Calculate the dydx for the graph
    dydx = finite_difference(x, y)

    # Calculate the running average for the found dydx
    dydx_average = running_average(dydx, n=n)

    start_idx = dydx_average.shape[0]//2

    for i in range(0, start_idx):
        if x[i] > 0.3:
            return i


def graph_alignment_2(x, y, tol=1e-8, n=11, stop=3600, num=14400, delta_x=0.1):

    num_cols = y.shape[1]

    # Arrays with the indices of the start and the plateau points
    indices_start = np.zeros(num_cols, dtype=np.int32)
    indices_plat = np.zeros(num_cols, dtype=np.int32)

    # Arrays with all the x values for start and the plateau points, and the time interval between those
    x_plat_values = np.zeros(num_cols, dtype=np.float64)
    x_start_values = np.zeros(num_cols, dtype=np.float64)
    x_start_to_plat = np.zeros(num_cols, dtype=np.float64)

    f = []  # List with interpolated functions
    for i in range(num_cols):
        # TODO at the ends boundary effects check if this is not a problem
        y[:, i] = running_average(y[:, i], n=n, mode="same")
        f_i = interp1d(x[:, i], y[:, i])

        x_i_interp1d = np.linspace(0, stop, num)
        y_i_interp1d = f_i(x_i_interp1d)

        # Calculate the start indices
        idx = find_start_index(
            x_i_interp1d, y_i_interp1d, tol=tol, n=n)
        indices_start[i] = idx
        x_start_values[i] = x[idx, i]

        # Calculate the plateau indices
        idx = find_platau_index(
            x_i_interp1d, y_i_interp1d, tol=tol, n=n)
        indices_plat[i] = idx
        x_plat_values[i] = x[idx, i]

        # Calculate the time between start and plateau
        x_start_to_plat = x_plat_values[i] - x_start_values[i]
        f.append(f_i)

    min_x_start_to_plat = np.amin(x_start_to_plat, axis=0)

    # Determine boundaries for the x domain
    x_boundaries = np.array(2, num_cols, dtype=np.float64)
    for i in range(num_cols):
        x_boundaries[1, i] = x_plat_values[i] - min_x_start_to_plat
        x_boundaries[2, i] = x_plat_values[i]

    x_shifted = []
    y_shifted = []

    for i in range(num_cols):
        num_x_i = (x_boundaries[2, i] - x_boundaries[1, i]) / delta_x
        x_i_domain = np.linspace(x_boundaries[1, i], x_boundaries[2, i], num_x)
        x_i_shifted = x_i_domain - x_plat_values[i]
        y_i_shifted = f[i](x_i_domain)

        # Appending shifted graphs to the lists
        x_shifted.append(x_i_shifted)
        y_shifted.append(y_i_shifted)

    return x_shifted, y_shifted


if __name__ == "__main__":
    pass
