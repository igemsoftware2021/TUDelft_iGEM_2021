from morris_method import morris_analysis
from models_parallelized import run_simulations_no_aptamer
from morris_method import morris_datawriter
import numpy as np

# Defining the properties of the Morris sensitivity analysis
trajectories = 15
num_levels = 4
num_parameters = 12

# Parameters
# (#) denotes the position in the parameters array
k_ts = 6.7*10**-5     # (0)  Transcription rate [uM/s]
k_tl = 7.2*10**-5     # (1)  Enzyme translation rate [uM/s]
k_mat = 0.5*10**-3    # (2)  Maturation rate of beta-galactosidase [1/s]
k_cat = 5.14*10**1    # (3)  Catalytic rate of beta-galactosidase [1/s]
k_s = 8.5*10**-3      # (4)  Michaelis constant of transcription [μM]
# (5)  Scaling factor for the transcription resources [1/s]
kc_s = 1.8*10**-4
k_l = 65.8*10**-3     # (6)  Michaelis constant of translation [μM]
k_tlr = 6*10**-5      # (7) Michaelis constant of translation resources [-]
k_m = 50              # (8) Michaelis constant of beta-galactosidase [μM]
deg_mrna = 1.3*10**-5  # (9) Degradation rate of mRNA [1/s]
deg_tlr = 7.5*10**-5  # (10) Degradation rate of translation resources [1/s]
k_on = 1*10**-2       # (11)  Association rate of vitamin and umRNA [1/µMs]
k_off = 1*10**-2      # (12)  Dissociation rate of vitamin and umRNA [1/s]
k_c = (1/60)/10            # (13)  Cleaving rate of umRNA [1/s]
parameters = np.array([k_ts, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                       k_tlr, k_m, deg_mrna, deg_tlr, k_on, k_off, k_c])  # Array containing above parameters

lower_range = 0.1
upper_range = 10


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
s_i = 1000  # Initial substrate concentration
initial_conditions = np.array([dna_i, s_i])

# Defining timescale of the model
t_tot = 1.5 * 3.600  # Total time [s]
dt = 0.01  # Timestep [s]


# Problem definition for prokaryotic system
no_aptamer_prokaryotic_problem = {
    'num_vars': num_parameters,
    'names': ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
              "k_tlr", "k_m", "deg_mrna", "deg_tlr", "dna_i"],
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
               [lower_range * dna_i, upper_range * dna_i]]  # (11) dna_i
}

# Eukaryotic problem definition (ignore this never used it)
no_aptamer_eukaryotic_problem = {
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
               [-1, 1],  # (10) deg_tlr
               [-1, 1]]  # (11) dna_i
}


# Doing Morris sensitivity analysis
(time, mu, mu_star, sigma, mu_star_conf_level) = morris_analysis(no_aptamer_prokaryotic_problem, trajectories,
                                                                 run_simulations_no_aptamer, constants, initial_conditions, dt=dt, t_tot=t_tot, num_levels=num_levels, seed=12)

# Saving the data
morris_datawriter(no_aptamer_prokaryotic_problem, "modelling\data\morris_no_aptamer",
                  "test", time, mu, mu_star, sigma, mu_star_conf_level)
