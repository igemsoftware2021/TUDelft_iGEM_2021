# Libraries to import
import numpy as np
import matplotlib.pyplot as plt

tbd = 1  # Placeholder for unkown parameters/concentrations so that they're easy to find

# Parameters
k_ts = tbd            # (0)  Transcription rate
k_on = 1              # (1)  Association rate of vitamin and umRNA [1/s]
k_off = 1             # (2)  Dissociation rate of vitamin and umRNA [1/s]
k_c = 1               # (3)  Cleaving rate of umRNA [s^-1]
k_tl = 0.05           # (4)  Enzyme translation rate [AA/s] !!!!!!!!!!!!!!!!
k_mat = 0.5*10**-3    # (5)  Maturation rate of beta-galactosidase [1/s]
k_cat = 4.28*10**-3   # (6)  Catalytic rate of beta-galactosidase [1/s]
k_s = 8.5*10**-3      # (7)  Michaelis constant of transcription [μM]
kc_s = 1.8*10**-4     # (8)  Scaling factor for the transcription resources [-]
k_l = 65.8*10**-3     # (9)  Michaelis constant of translation [μM]
k_tlr = 6*10**-6      # (10) Michaelis constant of translation resources [-]
k_m = 0.9             # (11) Michaelis constant of beta-galactosidase [μM]
deg_mrna = 1.4        # (12) Degradation rate of mRNA [1/s]
deg_tir = 7.5**10*-5  # (13) Degradation rate of translation resources [1/s]

parameters = np.array([k_ts, k_on, k_off, k_c, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                       k_tlr, k_m, deg_mrna, deg_tir])  # Array containing above parameters

t_tot = 3600  # total time [s]
dt = 0.1  # timestep [s]


# Initial concentrations
dna_i = tbd  # Initial concentration of the beta-galactosidase gene [μM]
vit_i = tbd   # Initial vitamin concentration
s_i = tbd   # Initial substrate concentration


def euk_model(parameters, dt, t_tot, dna_i, vit_i, s_i):
    n = np.ceil(t_tot/dt)  # number of timesteps of the simulation [-]

    # Arrays for the concentration of each species

    # Concentration of the beta-galactosidase gene [μM]
    dna = np.zeros(n, dtype=np.float64)

    # Concentration of uncleaved mRNA [μM]
    umrna = np.zeros(n, dtype=np.float64)

    # Concentration of cleaved mRNA [μM]
    cmrna = np.zeros(n, dtype=np.float64)

    # Concentration of degraged mRNA [μM]
    dmrna = np.zeros(n, dtype=np.float64)

    # Vitamine concentration [μM]
    vit = np.zeros(n, dtype=np.float64)

    # Concentration of the uncleaved mRNA - Vitamin complex [μM]
    umrna_vit = np.zeros(n, dtype=np.float64)

    # Prefactor that accounts for finite transcription resources and initiation rate [-]
    tsr = np.zeros(n, dtype=np.float64)

    # Prefactor that accounts for finite translation resources and initiation rate[-]
    tlr = np.zeros(n, dtype=np.float64)

    # Concentration of beta-galactosidase
    e = np.zeros(n, dtype=np.float64)
    s = np.zeros(n, dtype=np.float64)           #
    p = np.zeros(n, dtype=np.float64)           #

    # Fill in the initial concentrations
