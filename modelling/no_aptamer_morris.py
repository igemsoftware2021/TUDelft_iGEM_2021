from SALib.sample import morris as morris_sample
from no_aptamer_model import no_aptamer_model
from SALib.analyze import morris as morris_analyze
from SALib.test_functions import Ishigami
import numpy as np

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
n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-
time = np.linspace(0, t_tot, n)  # Array with all timepoints.

# Defining arrays for sensitivity indices. Each row contains the indices for one parameter over time.
mu = np.zeros((num_parameters, len(time)))  # The mean elementary effect
mu_star = np.zeros((num_parameters, len(time)))
sigma = np.zeros((num_parameters, len(time)))
mu_star_conf = np.zeros((num_parameters, len(time)))

# Generating input parameters for the model
model_input = morris_sample.sample(no_aptamer_problem, trajectories)

# to test only
# model_output = Ishigami.evaluate(model_input)
# indices = morris_analyze.analyze(no_aptamer_problem, model_input, model_output)
# print(indices)

# Number of different parameter sets to run the model with, so the number of simulations
simulations = np.shape(model_input)[0]
# An array for the output of the model (B/Y ratio)
model_output = np.zeros([len(time), simulations])

# Running the model for each set of parameters. The output is stored for each timepoint of each simulation
for ss in range(simulations):
    parameters = model_input[ss, :]
    model_output[:, ss] = no_aptamer_model(
        parameters, constants, dt, t_tot, dna_i, s_i)

# Running the Morris analysis at each timepoint (using tne output of all the different simulations)
for tt in range(len(time)):
    indices_dict = morris_analyze.analyze(no_aptamer_problem, model_input,
                                          model_output[tt, :], num_levels=4)
    mu[:, tt] = indices_dict['mu']
    mu_star[:, tt] = indices_dict['mu_star']
    sigma[:, tt] = indices_dict['sigma']
    mu_star_conf[:, tt] = indices_dict['mu_star_conf']
