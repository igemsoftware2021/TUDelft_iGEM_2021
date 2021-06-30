import numpy as np
from no_aptamer_model import no_aptamer_model
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from SALib.test_functions import Ishigami

no_aptamer_problem = {
    'num_vars': 11,
    'names': ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
              "k_tlr", "k_m", "deg_mrna", "deg_tlr"],
    'bounds': [[-1, 1],  # k_ts
               [-1, 1],  # k_tl
               [-1, 1],  # k_mat
               [-1, 1],  # k_cat
               [-1, 1],  # k_s
               [-1, 1],  # kc_s
               [-1, 1],  # k_l
               [-1, 1],  # k_tlr
               [-1, 1],  # k_m
               [-1, 1],  # deg_mrna
               [-1, 1], ]  # deg_tlr
}

trajectories = 10
levels = 4
dt = 1
t_tot = 100
dna_i = 1
s_i = 1

# Generate input parameters for the model
model_input = morris_sample.sample(no_aptamer_problem, trajectories)

# Determine the timepoints of the simulation
n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
time = np.linspace(0, t_tot, n)  # Array with all timepoints

# Simulations is the number of different parameter sets to run the model with, so the number of simulations
simulations = np.shape(model_input)[0]
model_output = np.zeros([len(time), simulations)

# Run the model for each set of parameters. The output is stored for each timepoint of each simulation
for ss in range(simulations):
    parameters = model_input[ss, :]
    model_output[:, ss] = no_aptamer_model(parameters, dt, t_tot, dna_i, s_i)


# put everything in morris analyzer,
# morris_analyze.analyze(no_aptamer_problem, model_input, model_output, num_levels=4)


# this would be morris with 10 trajectories and 4 levels and 11 parameters
