import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import time


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
    delta_x = x[1] - \
        x[0]  # TODO mogelijkheid van verschillende delta x binnen deze functie?
    # Calculate the forward difference for the first few points
    temp_forward = np.concatenate((
        np.array([0, 0], dtype=np.float64), y[:4]), axis=0)
    dydx_forward = np.convolve(
        temp_forward, np.array([-1, 4, -3, 0, 0], dtype=np.float64)/(2*delta_x), mode="valid")
    # Calculate the central difference for all the middle data points
    dydx_central = np.convolve(y, np.array(
        [-1, 8, 0, -8, 1], dtype=np.float64)/(12*delta_x), mode="valid")

    # Calculate the backward difference for the last few data points
    temp_backward = np.concatenate(
        (y[-4:], np.array([0, 0], dtype=np.float64)), axis=0)
    dydx_backward = np.convolve(temp_backward, np.array(
        [0, 0, 3, -4, 1], dtype=np.float64)/(2*delta_x), mode="valid")

    # Combine the computed finite differences and return them
    return np.concatenate((dydx_forward, dydx_central, dydx_backward), axis=0)


def find_platau_index(x, y, tol=1e-8, n=11):

    # Calculate the dydx for the graph
    dydx = finite_difference(x, y)

    # Calculate the running average for the found dydx
    dydx_average = running_average(dydx, n=n)

    start_idx = dydx_average.shape[0]//2

    for i in range(start_idx, dydx_average.shape[0]):
        if dydx_average[i] < tol:
            return i


def graph_alignment(x, y, tol=1e-8, n=11):
    """
    Function aligns all the graphs

    x: np.ndarray of shape (rows, columns), where every column is a different line
    """
    if len(x.shape) != 2 or len(y.shape) != 2:
        raise ValueError("Array x and array y should be two dimensional")

    num_cols = y.shape[1]

    indices = np.zeros(num_cols, dtype=np.int32)
    # Store all the x values for the plateau points
    x_plat_values = np.zeros(num_cols, dtype=np.float64)

    for i in range(num_cols):
        idx = find_platau_index(
            x[:, i], y[:, i], tol=tol, n=n)
        indices[i] = idx
        x_plat_values[i] = x[idx, i]

    idx_max_x_plat_value = np.argmax(x_plat_values, axis=0)
    max_x_plat_value = y[indices[idx_max_x_plat_value], idx_max_x_plat_value]

    print(idx_max_x_plat_value)

    x_shifted = x[:, idx_max_x_plat_value]
    y_shifted = np.zeros(y.shape, dtype=np.float64)

    # Find for every column what index is closest to the
    # found x platau value
    for i in range(num_cols):
        absolute_val_array = np.abs(y[:, i] - max_x_plat_value)
        idx = np.argmin(absolute_val_array)
        # The number of indices you need to shift
        idx_shift = idx_max_x_plat_value - idx
        if idx_shift == 0:
            y_shifted[:, i] = y[:, i]
        else:
            temp = y[:-idx_shift, i]
            y_shifted[idx_shift:, i] = y[:-idx_shift, i]

    return x_shifted, y_shifted


if __name__ == "__main__":
    pass
