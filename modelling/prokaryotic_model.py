# Libraries to import
import numpy as np
import matplotlib.pyplot as plt
from numba import njit


# @njit
def prok_model(parameters, dt, t_tot, dna_i, vit_i, s_i):
    # "Unpacking" the array with parameters into individual parameters
    k_ts = parameters[0]
    k_on = parameters[1]
    k_off = parameters[2]
    k_c = parameters[3]
    k_tl = parameters[4]
    k_mat = parameters[5]
    k_cat = parameters[6]
    k_s = parameters[7]
    kc_s = parameters[8]
    k_l = parameters[9]
    k_tlr = parameters[10]
    k_m = parameters[11]
    deg_mrna = parameters[12]
    deg_tlr = parameters[13]

    # Determine the timepoints of the simulation
    n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints

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
    # Concentration of monomeric subunits of beta-galactosidase
    e_mon = np.zeros(n, dtype=np.float64)
    # Concentration of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)

    # Plug in the initial concentrations
    dna[0] = dna_i
    vit[0] = vit_i
    s[0] = s_i
    tsr[0] = 1
    tlr[0] = 1

    # A loop with the differential equations
    for step in range(len(time)-1):
        # Differential of each species w.r.t time
        dna_dt = 0  # could remove this one, zero anyway
        umrna_dt = k_ts * tsr[step] * dna[step] / \
            (k_s + dna[step]) - k_c * umrna[step] - k_on * umrna[step] * \
            vit[step] + k_off * umrna_vit[step] - deg_mrna * umrna[step]
        cmrna_dt = k_c * umrna[step] - deg_mrna * cmrna[step]
        dmrna_dt = deg_mrna * (umrna[step] + cmrna[step])
        vit_dt = - k_on * vit[step] * umrna[step] + k_off * umrna_vit[step]
        umrna_vit_dt = - vit_dt
        tsr_dt = - kc_s * tsr[step] * dna[step] / (k_s + dna[step])
        tlr_dt = - deg_tlr * tlr[step] / (k_tlr + tlr[step])
        e_mon_dt = k_tl * tlr[step] * \
            cmrna[step] / (k_l + cmrna[step]) - \
            0.25 * k_mat * e_mon[step]
        e_dt = 0.25 * k_mat * e_mon[step]
        s_dt = - k_cat * e[step] * s[step] / (k_m + s[step])
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        dna[step + 1] = dna[step] + dna_dt * dt
        umrna[step + 1] = umrna[step] + umrna_dt * dt
        cmrna[step + 1] = cmrna[step] + cmrna_dt * dt
        dmrna[step + 1] = dmrna[step] + dmrna_dt * dt
        vit[step + 1] = vit[step] + vit_dt * dt
        umrna_vit[step + 1] = umrna_vit[step] + umrna_vit_dt * dt
        tsr[step + 1] = tsr[step] + tsr_dt * dt
        tlr[step + 1] = tlr[step] + tlr_dt * dt
        e_mon[step + 1] = e_mon[step] + e_mon_dt * dt
        e[step + 1] = e[step] + e_dt * dt
        s[step + 1] = s[step] + s_dt * dt
        p[step + 1] = p[step] + p_dt * dt
    return (p, time)


if __name__ == "__main__":
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
    dna_i = tbd  # Initial concentration of the beta-galactosidase gene [μM]
    vit_i = tbd   # Initial vitamin concentration
    s_i = tbd   # Initial substrate concentration

    # Running the model
    (data, time) = prok_model(parameters, dt, t_tot, dna_i, vit_i, s_i)
    (data_2, time_2) = prok_model_no_aptamer(parameters, dt, t_tot, dna_i, s_i)

    # Plotting to check if sth happened at all
    plt.plot(time, data)
    plt.show()
