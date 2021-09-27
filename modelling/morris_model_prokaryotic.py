import csv
import numpy as np
from morris_method import morris_analysis, morris_analysis_area
from morris_method import morris_datawriter
from models_area import model_prokaryotic_absorbance_area, model_prokaryotic_absorbance_area_parallel
from standard_values import standard_parameters_prokaryotic, standard_initial_conditions, standard_constants
import matplotlib.pyplot as plt

# Defining the properties of the Morris sensitivity analysis
trajectories = 30
num_levels = 4
num_parameters = 15

# Parameters
parameters = standard_parameters_prokaryotic()

# Parameter ranges
lower_range = 0.5
upper_range = 2


# Constants UNUSED ONLY LOOK AT CONCENTRATION CPR
# (#) denotes the position in the constants array
h = 0.020  # (0) Height of the paper [cm]
# (1) Exctinction coefficient of CPRG at a wavelength of 410 in 25% serum [1/(µM*cm)]
eps_cprg = 0.294
# (2) Exctinction coefficient of CPR at a wavelength of 580 in 25% serum [1/(µM*cm)]
eps_cpr = 0.539
# Array containing above constants
constants = np.array([h, eps_cprg, eps_cpr])

# Initial concentrations
dna_i = 5*10**-3  # Initial concentration of the beta-galactosidase gene [μM]
s_i = 1000  # Initial substrate concentration [μM]
vit_i = 14*10**-3  # Initial vitamin concentration [μM]
initial_conditions = np.array([dna_i, s_i, vit_i])

initial_conditions = standard_initial_conditions(
    dna_conc=5*10**-3, s_i=250, vit_conc=1)

# Defining timescale of the model
t_tot = 10800  # Total time [s]
dt = 0.01  # Timestep [s]

dna = 2*10**-3

# Problem definition for prokaryotic system
prokaryotic_problem = {
    'num_vars': num_parameters,
    'names': ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
              "k_tlr", "k_m", "deg_mrna", "deg_tlr", "k_on", "k_off", "k_c", "dna_i"],
    'bounds': [[lower_range * parameters[0], upper_range * parameters[0]],     # (0) k_ts
               [lower_range * parameters[1], upper_range * \
                   parameters[1]],     # (1) k_tl
               [lower_range * parameters[2], upper_range * \
                   parameters[2]],     # (2) k_mat
               [lower_range * parameters[3], upper_range * \
                   parameters[3]],     # (3) k_cat
               [lower_range * parameters[4], upper_range * \
                   parameters[4]],     # (4) k_s
               [lower_range * parameters[5], upper_range * \
                   parameters[5]],     # (5) kc_s
               [lower_range * parameters[6], upper_range * \
                   parameters[6]],     # (6) k_l
               [lower_range * parameters[7], upper_range * \
                   parameters[7]],     # (7) k_tlr
               [lower_range * parameters[8], upper_range * \
                   parameters[8]],     # (8) k_m
               [lower_range * parameters[9], upper_range * \
                   parameters[9]],     # (9) deg_mrna
               [lower_range * parameters[10], upper_range * \
                   parameters[10]],  # (10) deg_tlr
               [lower_range * parameters[11], upper_range * \
                   parameters[11]],  # (11) k_on
               [lower_range * parameters[12], upper_range * \
                   parameters[12]],  # (12) k_off
               [lower_range * parameters[13], upper_range * \
                   parameters[13]],  # ,  # (13) k_c
               [lower_range * dna_i, upper_range * dna_i]]  # (14) dna_i
    # TODO think about how to implement the changing dna concentration
    #    [lower_range * dna_i, upper_range * dna_i]]

}

# # Doing Morris sensitivity analysis
# (time, mu, mu_star, sigma, mu_star_conf_level) = morris_analysis(prokaryotic_problem, trajectories,
#                                                                  run_simulations_prokaryotic, constants, initial_conditions, dt=dt, t_tot=t_tot, num_levels=num_levels)

# (problem, trajectories, func, constants, initial_conditions, vit_conc1, vit_conc2, dt: int=0.01, t_tot: int=7200, num_levels: int=4, optimal_trajectories: int=None, local_optimization: bool=True, num_resamples: int=1000, conf_level: float=0.95, print_to_console: bool=False, seed: int=None):
# Doing Morris sensitivity analysis
(time, mu, mu_star, sigma, mu_star_conf_level) = morris_analysis_area(prokaryotic_problem, trajectories,
                                                                      model_prokaryotic_absorbance_area_parallel, constants=constants,
                                                                      dna_conc=5*10**-3, s_i=150, vit_conc1=0.05, vit_conc2=0.09,
                                                                      dt=dt, t_tot=t_tot, num_levels=num_levels)
# print(mu)
# plt.scatter(np.arange(0, mu.shape[0]), mu)

fig, ax = plt.subplots()
ax.bar(np.arange(0, mu.shape[0]), mu)
# plt.bar(np.arange(0, mu_star.shape[0]), mu_star)
# plt.hist(mu_star)
fig2, ax2 = plt.subplots()
ax2.bar(np.arange(0, sigma.shape[0]), sigma)
# plt.hist(sigma)
# print(mu_star_conf_level)
plt.show()

# Saving the data
path = "modelling\data\morris_prokaryotic"
tag = "test"
morris_datawriter(prokaryotic_problem, path,
                  tag, time, mu, mu_star, sigma, mu_star_conf_level)
# The names of the files are in the format:
# path/time_tag.csv
# path/mu_tag.csv
# etc.
