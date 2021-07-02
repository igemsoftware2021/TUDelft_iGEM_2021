from morris_method import morris_analysis
from models_parallelized import run_simulations_eukaryotic
import numpy as np
import csv

# Defining the properties of the Morris sensitivity analysis
trajectories = 10
num_levels = 4
num_parameters = 14

# Prokaryotic problem definition
prokaryotic_problem = {
    'num_vars': num_parameters,
    'names': ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
              "k_tlr", "k_m", "deg_mrna", "deg_tlr", "k_on", "k_off", "k_c"],
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
               [-1, 1],  # (10) deg_tlr
               [-1, 1],  # (11) k_on
               [-1, 1],  # (12) k_off
               [-1, 1], ]  # (13) k_c

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
v_i = 1  # Initial vitamin concentration
initial_conditions = np.array([dna_i, s_i])

# Defining timescale of the model
t_tot = 1  # Total time [s]
dt = 0.1  # Timestep [s]

# Doing Morris sensitivity analysis
(time, mu, mu_star, sigma, mu_star_conf_level) = morris_analysis(prokaryotic_problem, trajectories,
                                                                 run_simulations_eukaryotic, constants, initial_conditions, dt=dt, t_tot=t_tot, num_levels=num_levels)
