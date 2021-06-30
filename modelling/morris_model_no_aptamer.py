from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
import numpy as np
from numba import njit, prange

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
