
from models import model_no_aptamer
from models import model_eukaryotic
from models import model_prokaryotic
import numpy as np
from numba import njit, prange


@njit(parallel=True)
def run_simulations_no_aptamer(parameters, constants, initial_conditions, dt: int = 0.01, t_tot: int = 7200):
    """Function runs a simulations of the kinetics of the AptaVita system 
    without aptamers parallel using the Numba library.

    Parameters
    ----------
    parameters: numpy.array
        The Numpy array containing all the parameters for the model of dtype=float
    constants: numpy.array
        The Numpy array containing all the constants for the model of dtype=float
    initial_condition: numpy.array
        The Numpy array containing all the initial conditions for the model of dtype=float
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    -------
    model_output: numpy.ndarray
        The Numpy array contains all the model output values. Every row is a time point and
        every column contains the results of each unique simulation over all the timepoints.
        The resulting array has ceil(t_tot/dt) + 1 rows and parameters.shape[0] columns.
    """
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]
    model_output = np.zeros((n, num_simulations))
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        _, model_output[:, ii] = model_no_aptamer(
            parameters[ii, :], constants, initial_conditions, dt=dt, t_tot=t_tot)
    return model_output


@njit(parallel=True)
def run_simulations_prokaryotic(parameters, constants, initial_conditions, dt: int = 0.01, t_tot: int = 7200):
    """"Function runs a simulations of the kinetics of the prokaryotic
    AptaVita system parallel using the Numba library.

    Parameters
    ----------
    parameters: numpy.array
        The Numpy array containing all the parameters for the model of dtype=float
    constants: numpy.array
        The Numpy array containing all the constants for the model of dtype=float
    initial_condition: numpy.array
        The Numpy array containing all the initial conditions for the model of dtype=float
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    -------
    model_output: numpy.ndarray
        The Numpy array contains all the model output values. Every row is a time point and
        every column contains the results of each unique simulation over all the timepoints.
        The resulting array has ceil(t_tot/dt) + 1 rows and parameters.shape[0] columns.
    """
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]
    model_output = np.zeros((n, num_simulations))
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        _, model_output[:, ii] = model_prokaryotic(
            parameters[ii, :], constants, initial_conditions, dt=dt, t_tot=t_tot)
    return model_output


@njit(parallel=True)
def run_simulations_eukaryotic(parameters, constants, initial_conditions, dt: int = 0.01, t_tot: int = 7200):
    """Function runs a simulations of the kinetics of the prokaryotic
    AptaVita system parallel using the Numba library.

    Parameters
    ----------
    parameters: numpy.array
        The Numpy array containing all the parameters for the model of dtype=float
    constants: numpy.array
        The Numpy array containing all the constants for the model of dtype=float
    initial_condition: numpy.array
        The Numpy array containing all the initial conditions for the model of dtype=float
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    -------
    model_output: numpy.ndarray
        The Numpy array contains all the model output values. Every row is a time point and
        every column contains the results of each unique simulation over all the timepoints.
        The resulting array has ceil(t_tot/dt) + 1 rows and parameters.shape[0] columns.
    """
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]
    model_output = np.zeros((n, num_simulations))
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        _, model_output[:, ii] = model_eukaryotic(
            parameters[ii, :], constants, initial_conditions, dt=dt, t_tot=t_tot)
    return model_output
