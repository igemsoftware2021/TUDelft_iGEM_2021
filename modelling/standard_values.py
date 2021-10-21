import numpy as np


def standard_parameters_prokaryotic():
    # Parameters
    # (#) denotes the position in the parameters array
    k_ts = 6.7*10**-5     # (0)  Transcription rate [uM/s]
    k_tl = 7.0*10**-5     # (1)  Enzyme translation rate [uM/s]
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
    k_on = 2*10**-1  # (11)  Association rate of vitamin and umRNA [1/µMs]
    k_off = 1*10**-2      # (12)  Dissociation rate of vitamin and umRNA [1/s]
    k_c = 0.017         # (13)  (1/60) Cleaving rate of umRNA [1/s]
    parameters = np.array([k_ts, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                           k_tlr, k_m, deg_mrna, deg_tlr, k_on, k_off, k_c])  # Array containing above parameters

    return parameters


def standard_constants():
    # Constants
    # (#) denotes the position in the constants array
    h = 0.020  # (0) Height of the paper [cm]
    # (1) Exctinction coefficient of CPRG at a wavelength of 410 in 25% serum [1/(µM*cm)]
    eps_cprg = 0.294
    # (2) Exctinction coefficient of CPR at a wavelength of 574 nm [1/(µM*cm)] (Sigma-Aldrich)
    eps_cpr = 0.07228

    # Array containing above constants
    constants = np.array([h, eps_cprg, eps_cpr])
    return constants


def standard_initial_conditions():
    dna_conc = 3*10**-3     # DNA concentration [µM]
    s_i = 1000              # Substrate concentration [µM]
    vit_conc = 0.07         # Vitamin concentration [µM]
    return np.array([dna_conc, s_i, vit_conc])
