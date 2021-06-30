from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
import numpy as np
from numba import njit, prange


@njit(parallel=True)
def morris_run_simulations(func, parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
    """Function runs model functions """

    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]
    # The time steps
    model_output = np.zeros((n, num_simulations))
    for ii in prange(parameters.shape[0]):
        model_output[:, ii] = func(
            parameters[ii, :], constants, initial_conditions, dt=dt, t_tot=t_tot)
    return model_output


def morris_analysis(func, problem, trajectories, constants, initial_conditions, num_levels=4, dt=0.01, t_tot=7200):

    # Determening the timepoints of the simulation.
    n = int(np.ceil(t_tot/dt) + 1)  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints.

    # Defining arrays for sensitivity indices.
    # Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
    mu = np.zeros((n, num_parameters),
                  dtype=np.float32)  # The mean elementary effect
    mu_star = np.zeros((n, num_parameters), dtype=np.float32)
    sigma = np.zeros((n, num_parameters), dtype=np.float32)
    mu_star_conf_low = np.zeros((n, num_parameters), dtype=np.float32)
    mu_star_conf_high = np.zeros((n, num_parameters), dtype=np.float32)

    # Generating input parameters for the model
    # Each column is a parameter, one row contains all input parameters for one simulation
    model_input = morris_sample.sample(problem, trajectories)

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
        mu_star_conf_low[ii, :] = indices_dict['mu_star_conf'][0]
        mu_star_conf_high[ii, :] = indices_dict['mu_star_conf'][1]

    return time, mu, mu_star, sigma, mu_star_conf_low, mu_star_conf_high


# Defining the properties of the Morris sensitivity analysis
trajectories = 10
levels = 4
num_parameters = 11
no_aptamer_problem = {
    'num_vars': num_parameters,
    'names': ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
              "k_tlr", "k_m", "deg_mrna", "deg_tlr"],
    'bounds': [[-1, 1],  # (0) k_ts
               [-1, 1],  # (1) k_tl
               [-1, 1],  # (2) k_mat
               [-1, 1],  # (3) k_cat
               [-1, 1],  # (4) k_s
               [-1, 1],  # (5) kc_s
               [-1, 1],  # (6) k_l
               [-1, 1],  # (7) k_tlr
               [-1, 1],  # (8) k_m
               [-1, 1],  # (9) deg_mrna
               [-1, 1], ]  # (10) deg_tlr
}

# Constants
h = 8*10**-5  # (0) Height of the paper [cm]
eps_cprg = 1  # (1) Exctinction coefficient of CPRG at a wavelength of ???
eps_cpr = 1   # (2) Exctinction coefficient of CPR at a wavelength of ???
i0_cprg = 1   # (3) Blank measurement at a wavelength of ???
i0_cpr = 1    # (4) Blank measurement at a wavelength of ???
constants = np.array([h, eps_cprg, eps_cpr, i0_cprg, i0_cpr])

# Initial concentrations
dna_i = 1  # Initial concentration of the beta-galactosidase gene [μM]
s_i = 1   # Initial substrate concentration

# Defining timescale of the model
t_tot = 1  # Total time [s]
dt = 0.001  # Timestep [s]

# Determening the timepoints of the simulation.
n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
time = np.linspace(0, t_tot, n)  # Array with all timepoints.

# Defining arrays for sensitivity indices.
# Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
mu = np.zeros((time.shape[0], num_parameters),
              dtype=np.float32)  # The mean elementary effect
mu_star = np.zeros((time.shape[0], num_parameters), dtype=np.float32)
sigma = np.zeros((time.shape[0], num_parameters), dtype=np.float32)
mu_star_conf = np.zeros((time.shape[0], num_parameters), dtype=np.float32)

# Generating input parameters for the model
# Each column is a parameter, one row contains all input parameters for one simulation
model_input = morris_sample.sample(no_aptamer_problem, trajectories)

# Number of different parameter sets to run the model with, so the number of simulations
num_simulations = model_input.shape[0]
# An array for the output of the model (B/Y ratio)
model_output = np.zeros((time.shape[0], num_simulations))


# Running the model for each set of parameters. The output is stored for each timepoint of each simulation
for ii in range(simulations):
    parameters = model_input[ii, :]
    # Each column is a simulation, each row is the output of the simulation at one timestep
    model_output[:, ii] = no_aptamer_model(
        parameters, constants, dt, t_tot, dna_i, s_i)

# Running the Morris analysis at each timepoint (using tne output of all the different simulations)
for ii in range(time.shape[0]):
    indices_dict = morris_analyze.analyze(no_aptamer_problem, model_input,
                                          model_output[ii, :], num_levels=4)
    # Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
    mu[ii, :] = indices_dict['mu']
    mu_star[ii, :] = indices_dict['mu_star']
    sigma[ii, :] = indices_dict['sigma']
    mu_star_conf[ii, :] = indices_dict['mu_star_conf']
