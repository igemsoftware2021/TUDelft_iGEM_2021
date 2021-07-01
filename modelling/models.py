import numpy as np
from numba import njit

# @njit(cache=True, nogil=True)


def model_no_aptamer(parameters, constants, initial_conditions, dt=0.1, t_tot=7200):
    """Function does a simulation of the kinetics of the AptaVita system without aptamers.\n
    args:\n
    parameters: (np.ndarray)\n
    \tArray containing the model parameters (float): transcription rate, translation rate, maturation rate, 
        \tcatalytic rate of the enzyme, Michaelis constant transcription, scaling factor transcription resources, Michaelis constant translation, Michaelis constant translation resources, Michaelis constant enzymatic reaction, degradation rate mRNA, degredation rate transcription resources.\n
    constants: (np.ndarray)\n
    \tArray containing the model constants (float): height of the microfluidic paper, extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement CPRG, blank intensity measurement CPR.\n
    constants: (np.ndarray)\n
    \tArray containing the initial conditions: gene concentration, substrate concentration.\n
    dt: (int)\n
    \tThe time each timestep takes in seconds. (default 0.01)\n
    t_tot: (int)\n
    \tThe total time the model should run in seconds. (default 7200)\n
    \n
    return value:\n
    time, b_y: (tuple)\n
    \tTuple of which the first entry contains an array with all the timepoints of the simulation and the second entry the blue over yellow light intensity ratio at each timepoint.
    """
    # "Unpacking" the array with parameters into individual parameters
    k_ts = parameters[0]
    k_tl = parameters[1]
    k_mat = parameters[2]
    k_cat = parameters[3]
    k_s = parameters[4]
    kc_s = parameters[5]
    k_l = parameters[6]
    k_tlr = parameters[7]
    k_m = parameters[8]
    deg_mrna = parameters[9]
    deg_tlr = parameters[10]

    # "Unpacking" the array with constants into individual constants
    h = constants[0]
    eps_cprg = constants[1]
    eps_cpr = constants[2]
    i0_cprg = constants[3]
    i0_cpr = constants[4]

    # Determine the timepoints of the simulation
    n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints

    # Arrays for the concentration of each species
    # Concentration of the beta-galactosidase gene [μM]
    dna = np.zeros(n, dtype=np.float64)
    # Concentration of mRNA [μM]
    mrna = np.zeros(n, dtype=np.float64)
    # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
    dmrna = np.zeros(n, dtype=np.float64)
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

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna[0] = initial_conditions[0]
    s[0] = initial_conditions[1]
    tsr[0] = 1
    tlr[0] = 1

    # A loop with the differential equations
    for ii in range(time.shape[0]-1):
        # Differential of each species w.r.t time
        dna_dt = 0  # could remove this one, zero anyway
        mrna_dt = k_ts * tsr[ii] * dna[ii] / \
            (k_s + dna[ii]) - deg_mrna * mrna[ii]
        dmrna_dt = deg_mrna * mrna[ii]
        tsr_dt = - kc_s * tsr[ii] * dna[ii] / (k_s + dna[ii])
        tlr_dt = - deg_tlr * tlr[ii] / (k_tlr + tlr[ii])
        e_mon_dt = k_tl * tlr[ii] * \
            mrna[ii] / (k_l + mrna[ii]) - \
            k_mat * e_mon[ii]
        e_dt = k_mat * e_mon[ii]
        s_dt = - k_cat * e[ii] * s[ii] / (k_m + s[ii])
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        dna[ii + 1] = dna[ii] + dna_dt * dt
        mrna[ii + 1] = mrna[ii] + mrna_dt * dt
        dmrna[ii + 1] = dmrna[ii] + dmrna_dt * dt
        tsr[ii + 1] = tsr[ii] + tsr_dt * dt
        tlr[ii + 1] = tlr[ii] + tlr_dt * dt
        e_mon[ii + 1] = e_mon[ii] + e_mon_dt * dt
        e[ii + 1] = e[ii] + e_dt * dt
        s[ii + 1] = s[ii] + s_dt * dt
        p[ii + 1] = p[ii] + p_dt * dt

    # Calculating blue over yellow ratio
    blue = i0_cpr * 10 ** (eps_cpr * p * h)
    yellow = i0_cprg * 10 ** (eps_cprg * s * h)
    b_y = np.divide(blue, yellow)
    return time, b_y


# @njit(cache=True, nogil=True)
def model_prokaryotic(parameters, constants, initial_conditions, dt=0.1, t_tot=7200):
    """Function does a simulation of the kinetics of the prokaryotic AptaVita system.\n
    args:\n
    parameters: (np.ndarray)\n
    \tArray containing the model parameters (float): transcription rate, translation rate, maturation rate, 
        \tcatalytic rate of the enzyme, Michaelis constant transcription, scaling factor transcription resources, Michaelis constant translation, Michaelis constant translation resources, Michaelis constant enzymatic reaction, degradation rate mRNA, degredation rate transcription resources, association rate aptamer and vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.\n
    constants: (np.ndarray)\n
    \tArray containing the model constants (float): height of the microfluidic paper, extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement CPRG, blank intensity measurement CPR.\n
    constants: (np.ndarray)\n
    \tArray containing the initial conditions: gene concentration, substrate concentration.\n
    dt: (int)\n
    \tThe time each timestep takes in seconds. (default 0.01)\n
    t_tot: (int)\n
    \tThe total time the model should run in seconds. (default 7200)\n
    \n
    return value:\n
    (time, b_y): (tuple)\n
    \tTuple of which the first entry contains an array with all the timepoints of the simulation and the second entry the blue over yellow light intensity ratio at each timepoint.
    """
    # "Unpacking" the array with parameters into individual parameters
    k_ts = parameters[0]
    k_tl = parameters[1]
    k_mat = parameters[2]
    k_cat = parameters[3]
    k_s = parameters[4]
    kc_s = parameters[5]
    k_l = parameters[6]
    k_tlr = parameters[7]
    k_m = parameters[8]
    deg_mrna = parameters[9]
    deg_tlr = parameters[10]
    k_on = parameters[11]
    k_off = parameters[12]
    k_c = parameters[13]

    # "Unpacking" the array with constants into individual constants
    h = constants[0]
    eps_cprg = constants[1]
    eps_cpr = constants[2]
    i0_cprg = constants[3]
    i0_cpr = constants[4]

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

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna[0] = initial_conditions[0]
    s[0] = initial_conditions[1]
    vit[0] = initial_conditions[2]
    tsr[0] = 1
    tlr[0] = 1

    # A loop with the differential equations
    for step in range(time.shape[0]-1):
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

    # Calculating blue over yellow ratio
    blue = i0_cpr * 10 ** (eps_cpr * p * h)
    yellow = i0_cprg * 10 ** (eps_cprg * s * h)
    b_y = np.divide(blue, yellow)

    return time, b_y


# @njit(cache=True, nogil=True)
def model_eukaryotic(parameters, constants, initial_conditions, dt=0.1, t_tot=7200):
    """Function does a simulation of the kinetics of the eukaryotic AptaVita system.\n
    args:\n
    parameters: (np.ndarray)\n
    \tArray containing the model parameters (float): transcription rate, translation rate, maturation rate, 
        \tcatalytic rate of the enzyme, Michaelis constant transcription, scaling factor transcription resources, Michaelis constant translation, Michaelis constant translation resources, Michaelis constant enzymatic reaction, degradation rate mRNA, degredation rate transcription resources, association rate aptamer and vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.\n
    constants: (np.ndarray)\n
    \tArray containing the model constants (float): height of the microfluidic paper, extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement CPRG, blank intensity measurement CPR.\n
    constants: (np.ndarray)\n
    \tArray containing the initial conditions: gene concentration, substrate concentration.\n
    dt: (int)\n
    \tThe time each timestep takes in seconds. (default 0.01)\n
    t_tot: (int)\n
    \tThe total time the model should run in seconds. (default 7200)\n
    \n
    return value:\n
    (time, b_y): (tuple)\n
    \tTuple of which the first entry contains an array with all the timepoints of the simulation and the second entry the blue over yellow light intensity ratio at each timepoint.
    """
    # "Unpacking" the array with parameters into individual parameters
    k_ts = parameters[0]
    k_tl = parameters[1]
    k_mat = parameters[2]
    k_cat = parameters[3]
    k_s = parameters[4]
    kc_s = parameters[5]
    k_l = parameters[6]
    k_tlr = parameters[7]
    k_m = parameters[8]
    deg_mrna = parameters[9]
    deg_tlr = parameters[10]
    k_on = parameters[11]
    k_off = parameters[12]
    k_c = parameters[13]

    # "Unpacking" the array with constants into individual constants
    h = constants[0]
    eps_cprg = constants[1]
    eps_cpr = constants[2]
    i0_cprg = constants[3]
    i0_cpr = constants[4]

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

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna[0] = initial_conditions[0]
    s[0] = initial_conditions[1]
    vit[0] = initial_conditions[2]
    tsr[0] = 1
    tlr[0] = 1

    # A loop with the differential equations
    for ii in range(time.shape[0]-1):
        # Differential of each species w.r.t time
        dna_dt = 0  # could remove this one, zero anyway
        umrna_dt = k_ts * tsr[ii] * dna[ii] / \
            (k_s + dna[ii]) - k_c * umrna[ii] - k_on * umrna[ii] * \
            vit[ii] + k_off * umrna_vit[ii] - deg_mrna * umrna[ii]
        cmrna_dt = k_c * umrna[ii] - deg_mrna * cmrna[ii]
        dmrna_dt = deg_mrna * (umrna[ii] + cmrna[ii])
        vit_dt = - k_on * vit[ii] * umrna[ii] + k_off * umrna_vit[ii]
        umrna_vit_dt = - vit_dt
        tsr_dt = - kc_s * tsr[ii] * dna[ii] / (k_s + dna[ii])
        tlr_dt = - deg_tlr * tlr[ii] / (k_tlr + tlr[ii])
        e_mon_dt = k_tl * tlr[ii] * \
            (umrna_vit[ii] + umrna[ii]) / (k_l + (umrna_vit[ii] + umrna[ii])) - \
            k_mat * e_mon[ii]
        e_dt = k_mat * e_mon[ii]
        s_dt = - k_cat * e[ii] * s[ii] / (k_m + s[ii])
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        dna[ii + 1] = dna[ii] + dna_dt * dt
        umrna[ii + 1] = umrna[ii] + umrna_dt * dt
        cmrna[ii + 1] = cmrna[ii] + cmrna_dt * dt
        dmrna[ii + 1] = dmrna[ii] + dmrna_dt * dt
        vit[ii + 1] = vit[ii] + vit_dt * dt
        umrna_vit[ii + 1] = umrna_vit[ii] + umrna_vit_dt * dt
        tsr[ii + 1] = tsr[ii] + tsr_dt * dt
        tlr[ii + 1] = tlr[ii] + tlr_dt * dt
        e_mon[ii + 1] = e_mon[ii] + e_mon_dt * dt
        e[ii + 1] = e[ii] + e_dt * dt
        s[ii + 1] = s[ii] + s_dt * dt
        p[ii + 1] = p[ii] + p_dt * dt

    # Calculating blue over yellow ratio
    blue = i0_cpr * 10 ** (eps_cpr * p * h)
    yellow = i0_cprg * 10 ** (eps_cprg * s * h)
    b_y = np.divide(blue, yellow)

    return time, b_y
