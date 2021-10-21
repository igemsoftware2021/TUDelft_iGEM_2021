import numpy as np
from numpy.lib.function_base import diff
import scipy.signal
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from tqdm import tqdm


def weighted_running_average(x, y, n=11):  # n must be uneven
    num_x = x.shape[0]
    reference_index = int(np.floor(n/2))  # middle - 1
    y_avg = np.zeros((num_x - 2 * reference_index), dtype=np.float64)

    index_distance = np.arange(-reference_index,
                               reference_index+1, 1, dtype=np.int16)

    for i in tqdm(range(reference_index, x.shape[0] - reference_index)):
        weights = np.zeros(n, dtype=np.float64)
        for j in range(n):
            if j == reference_index:
                weights[j] = 1
            else:
                delta_x = abs(x[i] - x[i + index_distance[j]])
                weights[j] = abs(index_distance[j]) / delta_x
        total_weigth = np.sum(weights)
        y_avg[i-reference_index] = np.sum(np.multiply(
            y[i-reference_index:i+reference_index+1], weights)) / total_weigth
    x_avg = x[reference_index:(-reference_index)]
    return x_avg, y_avg


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

    start_idx = dydx.shape[0]//2
    for i in range(start_idx, dydx.shape[0]):
        if dydx[i] < tol:
           # print(i)
            return i + 1
# TODO so that code gives an error if no plateau
# TODO tolerance groter maken totdat je een plateau vindt


def find_start_index(x, y, tol=1e-8, n=11):
    # TODO change this whole function into sth that makes actual sense
    # Calculate the dydx for the graph
    dydx = finite_difference(x, y)

    start_idx = dydx.shape[0]//2

    for i in range(0, start_idx):
        if y[i] > 0.5:
            return i


def graph_alignment(x, y, tol=1e-8, n=11, num=14400, delta_x=0.01, kind="cubic"):

    num_cols = y.shape[1]

    # Arrays for storing the average values of y and the corresponding x datapoints
    num_boundary = int(np.floor(n/2))
    num_avg = x.shape[0] - 2 * num_boundary
    x_avg = np.zeros((num_avg, num_cols), dtype=np.float64)
    y_avg = np.zeros((num_avg, num_cols), dtype=np.float64)
    # Arrays for the indices of the start and the plateau points
    indices_start = np.zeros(num_cols, dtype=np.int32)
    indices_plat = np.zeros(num_cols, dtype=np.int32)

    # Arrays for all the x values for start and the plateau points, and the time interval between those
    x_plat_values = np.zeros(num_cols, dtype=np.float64)
    x_start_values = np.zeros(num_cols, dtype=np.float64)
    x_start_to_plat = np.zeros(num_cols, dtype=np.float64)

    f = []  # List with interpolated functions
    for i in range(num_cols):
        x_avg[:, i], y_avg[:, i] = weighted_running_average(
            x[:, i], y[:, i], n=n)
        f_i = interp1d(x_avg[:, i], y_avg[:, i], kind=kind)

        x_i_interp1d = np.linspace(x_avg[0, i], x_avg[-1, i], num)
        y_i_interp1d = f_i(x_i_interp1d)

        # Calculate the start indices
        idx = find_start_index(
            x_i_interp1d, y_i_interp1d, tol=tol, n=n)
        indices_start[i] = idx
        x_start_values[i] = x_i_interp1d[idx]
        # plt.figure(5)
        #plt.plot(x_i_interp1d, y_i_interp1d)
        # plt.show()
        # Calculate the plateau indices
        idx = find_platau_index(
            x_i_interp1d, y_i_interp1d, tol=tol, n=n)
        indices_plat[i] = idx
        x_plat_values[i] = x_i_interp1d[idx]

        # Calculate the time between start and plateau
        x_start_to_plat[i] = x_plat_values[i] - x_start_values[i]
        f.append(f_i)

    min_x_start_to_plat = np.amin(x_start_to_plat, axis=0)
    # print(x_start_values)
    # print(x_plat_values)
    # print(x_start_to_plat)
    # Determine boundaries for the x domain
    x_boundaries = np.zeros((2, num_cols), dtype=np.float64)
    for i in range(num_cols):
        x_boundaries[0, i] = x_plat_values[i] - min_x_start_to_plat
        x_boundaries[1, i] = x_plat_values[i]
    # print(x_boundaries)

    y_shifted = []

    # Shifting
    num_x = min_x_start_to_plat / delta_x
    num_x = np.int(np.floor(num_x))
    for i in range(num_cols):
        x_i_domain = np.linspace(
            x_boundaries[0, i], x_boundaries[1, i], num_x)
        if i == 0:
            x_shifted = x_i_domain - x_plat_values[i]
        y_i_shifted = f[i](x_i_domain)

        # Appending shifted graphs to the lists
        y_shifted.append(y_i_shifted)

    a_rmsd = np.zeros((num_cols, num_cols), dtype=np.float64)
    for i in range(num_cols):
        for j in range(num_cols):
            a_rmsd[i, j] = rmsd(x_shifted, y_shifted[i], y_shifted[j])
    return x_shifted, y_shifted, a_rmsd


def rmsd(x, y_1, y_2):
    # Calculating squared difference
    y_diff = y_2 - y_1
    y_sqrd_diff = y_diff * y_diff
    f = UnivariateSpline(x, y_sqrd_diff, k=1, s=0,
                         ext=2)  # Linear interpolation
    # Calculating mean of the squared difference
    y_def_int = f.integral(x[0], x[-1])
    x_interval = abs(x[-1] - x[0])
    y_mean_sqrd_diff = y_def_int / x_interval
    y_root_mean_sqrd_diff = y_mean_sqrd_diff ** 0.5  # Take root
    return y_root_mean_sqrd_diff


if __name__ == "__main__":
    pass
