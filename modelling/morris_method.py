from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
import numpy as np
from numba import njit, prange


@njit(parallel=True)
def morris_run_simulations(func, parameters, constants, initial_conditions, dt: int = 0.01, t_tot: int = 7200):
    """Function runs a function in parallel using the Numba library.\n
    \n
    Parameters\n
    ----------\n
    func: (function)\n
        The function that has as input args: parameters, constants, initial_conditions[, dt, t_tot]\n
    parameters: numpy.array\n
        The Numpy array containing all the parameters for the model of dtype=float\n
    constants: numpy.array\n
        The Numpy array containing all the constants for the model of dtype=float\n
    initial_condition: numpy.array\n
        The Numpy array containing all the initial conditions for the model of dtype=float\n
    dt: int\n
        The time each timestep takes in seconds. (default 0.01)\n
    t_tot: int\n
        The total time the model should run in seconds. (default 7200)\n
    \n
    Returns\n
    -------\n
    model_output: numpy.ndarray\n
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
    for ii in prange(parameters.shape[0]):
        _, model_output[ii, :] = func(
            parameters[ii, :], constants, initial_conditions, dt=dt, t_tot=t_tot)
    return model_output


def morris_analysis(problem, trajectories, func, constants, initial_conditions, dt: int = 0.01, t_tot: int = 7200, num_levels: int = 4, optimal_trajectories: int = None, local_optimization: bool = True, num_resamples: int = 1000, conf_level: float = 0.95, print_to_console: bool = False, seed: int = None):
    """Function does Morris analysis on a function/model.\n

    This is a wrapper function that creates model inputs required for Method of Morris. It then runs these model inputs through
    a model using a parallelized function. Then it analyses the model output and returns the Numpy arrays: time, mu, mu_star, sigma, mu_star_conf_level.
    \n
    This is a wrapper function that uses two functions from the SAlib library. Hence the description for the parameters:
    problem, trajectories, num_levels, optimal_trajectories, local_optimization, num_resamples, conf_level, print_to_console
    and seed where directly taken from https://salib.readthedocs.io/en/latest/api.html.

    \n
    Parameters\n
    ----------\n
    problem: dict\n
        The problem definition\n
    trajectories: int\n
        Number of trajectories to generate.\n
    func: (function)\n
        The function that has as input args: parameters, constants, initial_conditions[, dt, t_tot]\n
    parameters: numpy.array\n
        The Numpy array containing all the parameters for the model of dtype=float\n
    constants: numpy.array\n
        The Numpy array containing all the constants for the model of dtype=float\n
    initial_condition: numpy.array\n
        The Numpy array containing all the initial conditions for the model of dtype=float\n
    dt: int\n
        The time each timestep takes in seconds. (default 0.01)\n
    t_tot: int\n
        The total time the model should run in seconds. (default 7200)\n
    num_levels: int\n
        The number of grid levels (should be even) (default 4)\n
    optimal_trajectories: int\n
        The number of optimal trajectories to sample (between 2 and N) (default None)\n
    local_optimization: bool\n
        Flag whether to use local optimization according to Ruano et al. (2012) Speeds up the process
        tremendously for bigger N and num_levels. If set to False brute force method used, unless
        gurobipy is available (default True)\n
    num_resamples: int\n
        The number of resamples used to compute the confidence intervals (default 1000)\n
    conf_level: float\n
        The confidence interval level (default 0.95)
    print_to_console: bool\n
        Print results directly to console (default False)
    seed: int\n
        Seed to generate a random number (default None)\n
    \n
    Returns\n
    -------\n
    time: numpy.array


    model_output: (np.array (2D))\n
    Every row is a time point and the column contains the results of each unique simulation.

    """
    # First do some checks whether the problem variable contains all the needed keys for the function.
    if "num_vars" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'num_vars' (https://salib.readthedocs.io/en/latest/basics.html)")
    if "names" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'names' (https://salib.readthedocs.io/en/latest/basics.html)")
    if "bounds" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'bounds' (https://salib.readthedocs.io/en/latest/basics.html)")

    # Determening the timepoints of the simulation.
    n = int(np.ceil(t_tot/dt) + 1)  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints.

    # Retrieve number of parameters
    num_parameters = problem["num_vars"]

    # Defining arrays for sensitivity indices.
    # Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
    mu = np.zeros((n, num_parameters),
                  dtype=np.float32)  # The mean elementary effect
    mu_star = np.zeros((n, num_parameters), dtype=np.float32)
    sigma = np.zeros((n, num_parameters), dtype=np.float32)
    mu_star_conf_level = np.zeros((n, num_parameters), dtype=np.float32)

    # Generating input parameters for the model
    # Each column is a parameter, one row contains all input parameters for one simulation
    model_input = morris_sample.sample(
        problem, trajectories, num_levels=num_levels, optimal_trajectories=optimal_trajectories, local_optimization=local_optimization, seed=seed)

    model_output = morris_run_simulations(
        func, model_input, constants, initial_conditions, dt=dt, t_tot=t_tot)

    # Running the Morris analysis at each timepoint (using tne output of all the different simulations)
    for ii in range(n):
        indices_dict = morris_analyze.analyze(problem, model_input,
                                              model_output[ii, :], num_levels=num_levels)
        # Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
        mu[ii, :] = indices_dict['mu']
        mu_star[ii, :] = indices_dict['mu_star']
        sigma[ii, :] = indices_dict['sigma']
        mu_star_conf_level[ii, :] = np.array(
            indices_dict['mu_star_conf'], dtype=np.float32)

    return time, mu, mu_star, sigma, mu_star_conf_level
