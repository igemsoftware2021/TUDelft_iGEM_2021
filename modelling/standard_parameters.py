import numpy as np


def standard_parameters(dna_conc=5*10**-3, vit_conc=5):
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
    # (10) Degradation rate of translation resources [1/s]
    deg_tlr = 7.5*10**-5
    k_on = 1*10**-2       # (11)  Association rate of vitamin and umRNA [1/µMs]
    k_off = 1*10**-2      # (12)  Dissociation rate of vitamin and umRNA [1/s]
    k_c = (1/60)/10            # (13)  Cleaving rate of umRNA [1/s]
    parameters = np.array([k_ts, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                           k_tlr, k_m, deg_mrna, deg_tlr, k_on, k_off, k_c])  # Array containing above parameters

    # Constants
    # (#) denotes the position in the constants array
    h = 0.020  # (0) Height of the paper [cm]
    # (1) Exctinction coefficient of CPRG at a wavelength of 410 in 25% serum [1/(µM*cm)]
    eps_cprg = 0.294
    # (2) Exctinction coefficient of CPR at a wavelength of 580 in 25% serum [1/(µM*cm)]
    eps_cpr = 0.539
    # I_null = ... Has to be filled in after the experiments.

    # Array containing above constants
    constants = np.array([h, eps_cprg, eps_cpr])

    # Initial concentration of the beta-galactosidase gene [μM]
    # dna_i = 5*10**-3
    s_i = 150         # Initial substrate concentration [μM]

    initial_conditions = np.array([dna_conc, s_i, vit_conc])
    return parameters, constants, initial_conditions
