# Libraries to import
from models import model_prokaryotic
import numpy as np
import matplotlib.pyplot as plt
import csv

tbd = 1  # Placeholder for unkown parameters/concentrations so that they're easy to find

# Parameters
# (#) denotes the position in the parameters array
k_ts = tbd            # (0)  Transcription rate
k_on = 1              # (1)  Association rate of vitamin and umRNA [1/s]
k_off = 1             # (2)  Dissociation rate of vitamin and umRNA [1/s]
k_c = 1               # (3)  Cleaving rate of umRNA [s^-1]
k_tl = 3*10**-5       # (4)  Enzyme translation rate [1/s]
k_mat = 0.5*10**-3    # (5)  Maturation rate of beta-galactosidase [1/s]
k_cat = 4.28*10**3   # (6)  Catalytic rate of beta-galactosidase [1/s]
k_s = 8.5*10**-3      # (7)  Michaelis constant of transcription [μM]
# (8)  Scaling factor for the transcription resources [-]
kc_s = 1.8*10**-4
k_l = 65.8*10**-3     # (9)  Michaelis constant of translation [μM]
# (10) Michaelis constant of translation resources [-]
k_tlr = 6*10**-6
k_m = 0.9             # (11) Michaelis constant of beta-galactosidase [μM]
deg_mrna = 1.3*10**-5  # (12) Degradation rate of mRNA [1/s]
# (13) Degradation rate of translation resources [1/s]
deg_tlr = 7.5*10**-5
parameters = np.array([k_ts, k_on, k_off, k_c, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                       k_tlr, k_m, deg_mrna, deg_tlr])  # Array containing above parameters

t_tot = 3600  # total time [s]
dt = 0.01  # timestep [s]

# Initial concentrations
dna_i = 5*10**-3  # Initial concentration of the beta-galactosidase gene [μM]
s_i = tbd   # Initial substrate concentration
vit_i = tbd   # Initial vitamin concentration
# Array containing above constants
initial_conditions = np.array([dna_i, s_i, vit_i])

# Constants
# (#) denotes the position in the constants array
h = 8*10**-5 * 1  # (0) Height of the paper [cm]
eps_cprg = tbd  # (1) Exctinction coefficient of CPRG at a wavelength of ???
eps_cpr = tbd   # (2) Exctinction coefficient of CPR at a wavelength of ???
# Array containing above constants
constants = np.array([h, eps_cprg, eps_cpr])

# Running the model
(time, data) = model_prokaryotic(parameters,
                                 constants, initial_conditions, dt=dt, t_tot=t_tot)

# Plotting
plt.plot(time, data)
plt.show()
