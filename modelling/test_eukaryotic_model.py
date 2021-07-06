import numpy as np
from numba import njit, prange

# https://stackoverflow.com/questions/16924688/using-timeit-module-in-a-function-with-arguments

# @njit


def euk_model(parameters, dt, t_tot, dna_i, vit_i, s_i):
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
    h = parameters[14]
    eps_cprg = parameters[15]
    eps_cpr = parameters[16]
    i0_cprg = parameters[17]
    i0_cpr = parameters[18]

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
    # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
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
    # Concentration of active sites of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow intensity ratio
    b_y = np.zeros(n, dtype=np.float64)

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
            (umrna_vit[step] + umrna[step]) / (k_l + (umrna_vit[step] + umrna[step])) - \
            k_mat * e_mon[step]
        e_dt = k_mat * e_mon[step]
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

    return time, p


@njit
def njit_euk_model(parameters, dt, t_tot, dna_i, vit_i, s_i):
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
    h = parameters[14]
    eps_cprg = parameters[15]
    eps_cpr = parameters[16]
    i0_cprg = parameters[17]
    i0_cpr = parameters[18]

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
    # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
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
    # Concentration of active sites of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow intensity ratio
    b_y = np.zeros(n, dtype=np.float64)

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
            (umrna_vit[step] + umrna[step]) / (k_l + (umrna_vit[step] + umrna[step])) - \
            k_mat * e_mon[step]
        e_dt = k_mat * e_mon[step]
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

    return time, p


@njit
def njit_euk_model_one(parameters, dt, t_tot, dna_i, vit_i, s_i):
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
    h = parameters[14]
    eps_cprg = parameters[15]
    eps_cpr = parameters[16]
    i0_cprg = parameters[17]
    i0_cpr = parameters[18]

    # Determine the timepoints of the simulation
    n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints

    results = np.zeros((13, n), dtype=np.float64)

    # Computing the concentration of each species for the next step
    # dna: results[0, ii]
    # umrna: results[1, ii]
    # cmrna: results[2, ii]
    # dmrna: results[3, ii]
    # vit: results[4, ii]
    # umrna_vit: results[5, ii]
    # tsr: results[6, ii]
    # tlr: results[7, ii]
    # e_mon: results[8, ii]
    # e: results[9, ii]
    # s: results[10, ii]
    # p: results[11, ii]

    # results[0, ii]
    # # Arrays for the concentration of each species
    # # Concentration of the beta-galactosidase gene [μM]
    # dna = np.zeros(n, dtype=np.float64)
    # # Concentration of uncleaved mRNA [μM]
    # umrna = np.zeros(n, dtype=np.float64)
    # # Concentration of cleaved mRNA [μM]
    # cmrna = np.zeros(n, dtype=np.float64)
    # # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
    # dmrna = np.zeros(n, dtype=np.float64)
    # # Vitamine concentration [μM]
    # vit = np.zeros(n, dtype=np.float64)
    # # Concentration of the uncleaved mRNA - Vitamin complex [μM]
    # umrna_vit = np.zeros(n, dtype=np.float64)
    # # Prefactor that accounts for finite transcription resources and initiation rate [-]
    # tsr = np.zeros(n, dtype=np.float64)
    # # Prefactor that accounts for finite translation resources and initiation rate[-]
    # tlr = np.zeros(n, dtype=np.float64)
    # # Concentration of monomeric subunits of beta-galactosidase
    # e_mon = np.zeros(n, dtype=np.float64)
    # # Concentration of active sites of beta-galactosidase (enzyme)
    # e = np.zeros(n, dtype=np.float64)
    # # Concentration of CPRG (substrate)
    # s = np.zeros(n, dtype=np.float64)
    # # Concentration of CPR (product)
    # p = np.zeros(n, dtype=np.float64)
    # # Blue over yellow intensity ratio
    # b_y = np.zeros(n, dtype=np.float64)

    # Plug in the initial concentrations
    results[0, 0] = dna_i
    results[4, 0] = vit_i
    results[10, 0] = s_i
    results[6, 0] = 1
    results[7, 0] = 1

    # A loop with the differential equations
    for ii in range(len(time)-1):
        # Differential of each species w.r.t time
        dna_dt = 0  # could remove this one, zero anyway
        umrna_dt = k_ts * results[6, ii] * results[0, ii+1] / \
            (k_s + results[0, ii+1]) - k_c * results[1, ii] - k_on * results[1, ii] * \
            results[4, ii] + k_off * results[5, ii] - deg_mrna * results[1, ii]
        cmrna_dt = k_c * results[1, ii] - deg_mrna * results[2, ii]
        dmrna_dt = deg_mrna * (results[1, ii] + results[2, ii])
        vit_dt = - k_on * results[4, ii] * \
            results[1, ii] + k_off * results[5, ii]
        umrna_vit_dt = - vit_dt
        tsr_dt = - kc_s * results[6, ii] * \
            results[0, ii] / (k_s + results[0, ii])
        tlr_dt = - deg_tlr * results[7, ii] / (k_tlr + results[7, ii])
        e_mon_dt = k_tl * results[7, ii] * \
            (results[5, ii] + results[1, ii]) / (k_l + (results[5, ii] + results[1, ii])) - \
            k_mat * results[8, ii]
        e_dt = k_mat * results[8, ii]
        s_dt = - k_cat * results[9, ii] * \
            results[10, ii] / (k_m + results[10, ii])
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        results[0, ii+1] = results[0, ii] + dna_dt * dt
        results[1, ii+1] = results[1, ii] + umrna_dt * dt
        results[2, ii+1] = results[2, ii] + cmrna_dt * dt
        results[3, ii+1] = results[3, ii] + dmrna_dt * dt
        results[4, ii+1] = results[4, ii] + vit_dt * dt
        results[5, ii+1] = results[5, ii] + umrna_vit_dt * dt
        results[6, ii+1] = results[6, ii] + tsr_dt * dt
        results[7, ii+1] = results[7, ii] + tlr_dt * dt
        results[8, ii+1] = results[8, ii] + e_mon_dt * dt
        results[9, ii+1] = results[9, ii] + e_dt * dt
        results[10, ii+1] = results[10, ii] + s_dt * dt
        results[11, ii+1] = results[11, ii] + p_dt * dt

    return time, results[-1, :]


@njit(cache=True)
def njit_euk_model_cache(parameters, dt, t_tot, dna_i, vit_i, s_i):
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
    h = parameters[14]
    eps_cprg = parameters[15]
    eps_cpr = parameters[16]
    i0_cprg = parameters[17]
    i0_cpr = parameters[18]

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
    # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
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
    # Concentration of active sites of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow intensity ratio
    b_y = np.zeros(n, dtype=np.float64)

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
            (umrna_vit[step] + umrna[step]) / (k_l + (umrna_vit[step] + umrna[step])) - \
            k_mat * e_mon[step]
        e_dt = k_mat * e_mon[step]
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

    return time, p


@njit(nogil=True)
def njit_euk_model_nogil(parameters, dna_i, vit_i, s_i, dt=0.1, t_tot=3600):
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
    h = parameters[14]
    eps_cprg = parameters[15]
    eps_cpr = parameters[16]
    i0_cprg = parameters[17]
    i0_cpr = parameters[18]

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
    # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
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
    # Concentration of active sites of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow intensity ratio
    b_y = np.zeros(n, dtype=np.float64)

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
            (umrna_vit[step] + umrna[step]) / (k_l + (umrna_vit[step] + umrna[step])) - \
            k_mat * e_mon[step]
        e_dt = k_mat * e_mon[step]
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

    return time, p


@njit()
def run_non_parallel_nogil(parameters, dna_i, vit_i, s_i, dt=0.01, t_tot=3600):
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]
    model_output = np.zeros((n, num_simulations))
    # Every column is a unique simulation
    for ii in range(num_simulations):
        _, model_output[:, ii] = njit_euk_model_nogil(
            parameters[ii, :], dna_i, vit_i, s_i, dt=dt, t_tot=t_tot)
    return model_output


@njit(parallel=True)
def run_parallel_nogil(parameters, dna_i, vit_i, s_i, dt=0.01, t_tot=3600):
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]
    model_output = np.zeros((n, num_simulations))
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        _, model_output[:, ii] = njit_euk_model_nogil(
            parameters[ii, :], dna_i, vit_i, s_i, dt=dt, t_tot=t_tot)
    return model_output


if __name__ == "__main__":
    import timeit
    tbd = 1  # Placeholder for unkown parameters/concentrations so that they're easy to find

    # Parameters
    # (#) denotes the position in the parameters array
    k_ts = tbd            # (0)  Transcription rate
    k_on = 1              # (1)  Association rate of vitamin and umRNA [1/s]
    k_off = 1             # (2)  Dissociation rate of vitamin and umRNA [1/s]
    k_c = 1               # (3)  Cleaving rate of umRNA [s^-1]
    k_tl = 5*10**-5       # (4)  Enzyme translation rate [1/s]
    k_mat = 0.5*10**-3    # (5)  Maturation rate of beta-galactosidase [1/s]
    k_cat = 4.28*10**3    # (6)  Catalytic rate of beta-galactosidase [1/s]
    k_s = 8.5*10**-3      # (7)  Michaelis constant of transcription [μM]
    # (8)  Scaling factor for the transcription resources [-]
    kc_s = 1.8*10**-4
    k_l = 65.8*10**-3     # (9)  Michaelis constant of translation [μM]
    # (10) Michaelis constant of translation resources [-]
    k_tlr = 6*10**-6
    k_m = 0.9             # (11) Michaelis constant of beta-galactosidase [μM]
    deg_mrna = 1.4*10**-3  # (12) Degradation rate of mRNA [1/s]
    # (13) Degradation rate of translation resources [1/s]
    deg_tlr = 7.5*10**-5
    h = 8*10**-5  # (14) Height of the paper [cm]
    eps_cprg = 1  # (15) Exctinction coefficient of CPRG at a wavelength of ???
    eps_cpr = 1   # (16) Exctinction coefficient of CPR at a wavelength of ???
    i0_cprg = 1   # (17) Blank measurement at a wavelength of ???
    i0_cpr = 1    # (18) Blank measurement at a wavelength of ???
    parameters = np.array([k_ts, k_on, k_off, k_c, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                           k_tlr, k_m, deg_mrna, deg_tlr, h, eps_cprg, eps_cpr, i0_cprg, i0_cpr])  # Array containing above parameters

    n_dim = 500
    parameter_matrix = np.zeros((n_dim, parameters.shape[0]), dtype=np.float32)
    for j in range(n_dim):
        parameter_matrix[j, :] = parameters

    t_tot = 3600  # total time [s]
    dt = 0.01  # timestep [s]

    # Initial concentrations
    dna_i = tbd  # Initial concentration of the beta-galactosidase gene [μM]
    vit_i = tbd   # Initial vitamin concentration
    s_i = tbd   # Initial substrate concentration

    # print(timeit.repeat("euk_model(parameters, dt, t_tot, dna_i, vit_i, s_i)",
    #       setup="from __main__ import euk_model, parameters, dt, t_tot, dna_i, vit_i, s_i", number=1, repeat=1))

    # print(timeit.repeat("njit_euk_model(parameters, dt, t_tot, dna_i, vit_i, s_i)",
    #       setup="from __main__ import njit_euk_model, parameters, dt, t_tot, dna_i, vit_i, s_i", number=100, repeat=3))

    # print(timeit.repeat("njit_euk_model_one(parameters, dt, t_tot, dna_i, vit_i, s_i)",
    #       setup="from __main__ import njit_euk_model_one, parameters, dt, t_tot, dna_i, vit_i, s_i", number=100, repeat=3))

    # print(timeit.repeat("njit_euk_model_cache(parameters, dt, t_tot, dna_i, vit_i, s_i)",
    #       setup="from __main__ import njit_euk_model_cache, parameters, dt, t_tot, dna_i, vit_i, s_i", number=100, repeat=3))

    # print(timeit.repeat("njit_euk_model_nogil(parameters, dt, t_tot, dna_i, vit_i, s_i)",
    #       setup="from __main__ import njit_euk_model_nogil, parameters, dt, t_tot, dna_i, vit_i, s_i", number=100, repeat=3))

    print(timeit.repeat("run_non_parallel_nogil(parameter_matrix, dna_i, vit_i, s_i, dt=dt, t_tot=t_tot)",
          setup="from __main__ import run_non_parallel_nogil, parameter_matrix, dt, t_tot, dna_i, vit_i, s_i", number=1, repeat=2))

    print(timeit.repeat("run_parallel_nogil(parameter_matrix, dna_i, vit_i, s_i, dt=dt, t_tot=t_tot)",
          setup="from __main__ import run_parallel_nogil, parameter_matrix, dt, t_tot, dna_i, vit_i, s_i", number=1, repeat=2))
