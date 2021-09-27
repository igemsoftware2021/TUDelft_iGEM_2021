import numpy as np
from numba import njit
import csv
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions


@njit(cache=True, nogil=True)
def model_yfp_expression(parameters, dna_i, dt=0.01, t_tot=7200):
    """Function does a simulation of YFP expression.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity 
        measurement CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate concentration.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
    """
    # "Unpacking" the array with parameters into individual parameters
    k_ts = parameters[0]
    k_tl = parameters[1]
    k_mat = parameters[2]
    k_s = parameters[3]
    kc_s = parameters[4]
    k_l = parameters[5]
    k_tlr = parameters[6]
    deg_mrna = parameters[7]
    deg_tlr = parameters[8]

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
    yfp = np.zeros(n, dtype=np.float64)
    # Concentration of beta-galactosidase (enzyme)
    yfp_mat = np.zeros(n, dtype=np.float64)

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna[0] = dna_i
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
        yfp_dt = k_tl * tlr[ii] * \
            mrna[ii] / (k_l + mrna[ii]) - \
            0.25 * k_mat * yfp[ii]
        yfp_mat_dt = 0.25 * k_mat * yfp[ii]

        # Computing the concentration of each species for the next step
        dna[ii + 1] = dna[ii] + dna_dt * dt
        mrna[ii + 1] = mrna[ii] + mrna_dt * dt
        dmrna[ii + 1] = dmrna[ii] + dmrna_dt * dt
        tsr[ii + 1] = tsr[ii] + tsr_dt * dt
        tlr[ii + 1] = tlr[ii] + tlr_dt * dt
        yfp[ii + 1] = yfp[ii] + yfp_dt * dt
        yfp_mat[ii + 1] = yfp_mat[ii] + yfp_mat_dt * dt
    return time, yfp_mat


@njit(cache=True, nogil=True)
def model_no_aptamer(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
    """Function does a simulation of the kinetics of the AptaVita system without aptamers.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity 
        measurement CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate concentration.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
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
    # Concentration of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow absorbance ratio
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
            0.25 * k_mat * e_mon[ii]
        e_dt = 0.25 * k_mat * e_mon[ii]
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
    blue = eps_cpr * p * h
    yellow = eps_cprg * s * h
    b_y = np.divide(blue, yellow)
    return time, p


@njit(cache=True, nogil=True)
def model_prokaryotic(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
    """Function does a simulation of the kinetics of the prokaryotic AptaVita system.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources association rate aptamer and 
        vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement 
        CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate 
        concentration, vitamin concentration.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
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

    cmrna_dt_array = np.zeros(n, dtype=np.float64)

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
    # Concentration of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow absorbance ratio
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
        cmrna_dt_array[step + 1] = cmrna_dt * 1/dt

    # Calculating blue over yellow ratio
    blue = eps_cpr * p * h
    yellow = eps_cprg * s * h
    b_y = np.divide(blue, yellow)
    return time, p


@njit(cache=True, nogil=True)
def model_prokaryotic_all(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
    """Function does a simulation of the kinetics of the prokaryotic AptaVita system.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources association rate aptamer and 
        vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement 
        CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate 
        concentration, vitamin concentration.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
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

    cmrna_dt_array = np.zeros(n, dtype=np.float64)

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
    # Concentration of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow absorbance ratio
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
        cmrna_dt_array[step + 1] = cmrna_dt * 1/dt

    # Calculating blue over yellow ratio
    blue = eps_cpr * p * h
    yellow = eps_cprg * s * h
    b_y = np.divide(blue, yellow)

    return time, dna, umrna, cmrna, dmrna, vit, umrna_vit, tsr, tlr, e_mon, e, s, p


@njit(cache=True, nogil=True)
def model_prokaryotic_hp(parameters, constants, initial_conditions, dt=0.01, t_tot=7200, dt_store=0.01):
    """Function does a high performance simulation of the kinetics of the prokaryotic AptaVita system.
    The simulation is run at dt timesteps, but the data is stored at dt_store timestep, since the size
    of the arrays will otherwise be to large.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources association rate aptamer and 
        vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement 
        CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate 
        concentration, vitamin concentration.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
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

    if dt >= dt_store:
        dt_store = dt

    # Determine the timepoints of the simulation
    # Number of timesteps to store of the simulation [-]
    store_steps = int(np.ceil(dt_store/dt))
    n = int(np.ceil((t_tot/dt))/store_steps) + 1
    time = np.linspace(0, t_tot, n)  # Array with all timepoints

    # Arrays for the concentration of each species
    # Concentration of the beta-galactosidase gene [μM]
    dna = 0
    dna_store = np.zeros(n, dtype=np.float64)
    # Concentration of uncleaved mRNA [μM]
    umrna = 0
    umrna_store = np.zeros(n, dtype=np.float64)
    # Concentration of cleaved mRNA [μM]
    cmrna = 0
    cmrna_store = np.zeros(n, dtype=np.float64)

    # Concentration of degraded mRNA [μM] (this does not denote a physical concentraion within the system, but merely tracks the concentration of mRNA that has been degraded)
    dmrna = 0
    dmrna_store = np.zeros(n, dtype=np.float64)
    # Vitamine concentration [μM]
    vit = 0
    vit_store = np.zeros(n, dtype=np.float64)
    # Concentration of the uncleaved mRNA - Vitamin complex [μM]
    umrna_vit = 0
    umrna_vit_store = np.zeros(n, dtype=np.float64)
    # Prefactor that accounts for finite transcription resources and initiation rate [-]
    tsr = 0
    tsr_store = np.zeros(n, dtype=np.float64)
    # Prefactor that accounts for finite translation resources and initiation rate[-]
    tlr = 0
    tlr_store = np.zeros(n, dtype=np.float64)
    # Concentration of monomeric subunits of beta-galactosidase
    e_mon = 0
    e_mon_store = np.zeros(n, dtype=np.float64)
    # Concentration of beta-galactosidase (enzyme)
    e = 0
    e_store = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = 0
    s_store = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = 0
    p_store = np.zeros(n, dtype=np.float64)
    # Blue over yellow absorbance ratio
    b_y_store = np.zeros(n, dtype=np.float64)

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna_store[0] = initial_conditions[0]
    s_store[0] = initial_conditions[1]
    vit_store[0] = initial_conditions[2]
    tsr_store[0] = 1
    tlr_store[0] = 1

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna = initial_conditions[0]
    s = initial_conditions[1]
    vit = initial_conditions[2]
    tsr = 1
    tlr = 1

    # A loop with the differential equations
    for step in range(int(np.ceil(t_tot/dt))):
        # Differential of each species w.r.t time
        dna_dt = 0  # could remove this one, zero anyway
        umrna_dt = k_ts * tsr * dna / \
            (k_s + dna) - k_c * umrna - k_on * umrna * \
            vit + k_off * umrna_vit - deg_mrna * umrna
        cmrna_dt = k_c * umrna - deg_mrna * cmrna
        dmrna_dt = deg_mrna * (umrna + cmrna)
        vit_dt = - k_on * vit * umrna + k_off * umrna_vit
        umrna_vit_dt = - vit_dt
        tsr_dt = - kc_s * tsr * dna / (k_s + dna)
        tlr_dt = - deg_tlr * tlr / (k_tlr + tlr)
        e_mon_dt = k_tl * tlr * \
            cmrna / (k_l + cmrna) - \
            0.25 * k_mat * e_mon
        e_dt = 0.25 * k_mat * e_mon
        s_dt = - k_cat * e * s / (k_m + s)
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        dna = dna + dna_dt * dt
        umrna = umrna + umrna_dt * dt
        cmrna = cmrna + cmrna_dt * dt
        dmrna = dmrna + dmrna_dt * dt
        vit = vit + vit_dt * dt
        umrna_vit = umrna_vit + umrna_vit_dt * dt
        tsr = tsr + tsr_dt * dt
        tlr = tlr + tlr_dt * dt
        e_mon = e_mon + e_mon_dt * dt
        e = e + e_dt * dt
        s = s + s_dt * dt
        p = p + p_dt * dt

        if step % store_steps == 0:
            index = step//store_steps
            dna_store[index] = dna
            umrna_store[index] = umrna
            cmrna_store[index] = cmrna
            dmrna_store[index] = dmrna
            vit_store[index] = vit
            umrna_vit_store[index] = umrna_vit
            tsr_store[index] = tsr
            tlr_store[index] = tlr
            e_mon_store[index] = e_mon
            e_store[index] = e
            s_store[index] = s
            p_store[index] = p

    # Calculating blue over yellow ratio
    blue = eps_cpr * p * h
    yellow = eps_cprg * s * h
    b_y = np.divide(blue, yellow)
    return time, p_store


@njit(cache=True, nogil=True)
def model_eukaryotic(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
    """Function does a simulation of the kinetics of the eukaryotic AptaVita system.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources association rate aptamer and 
        vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement 
        CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate 
        concentration, vitamin concentration.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
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
    # Concentration of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow absorbance ratio
    b_y = np.zeros(n, dtype=np.float64)

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
            0.25 * k_mat * e_mon[ii]
        e_dt = 0.25 * k_mat * e_mon[ii]
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
    blue = eps_cpr * p * h
    yellow = eps_cprg * s * h
    b_y = np.divide(blue, yellow)
    return time, b_y


@njit(cache=True, nogil=True)
def model_prokaryotic_readout(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
    """Function does a simulation of the kinetics of the prokaryotic AptaVita system for different vitamin concentrations.
    There is some timedelay in the readout to visualize differences in viscosity of blood plasma and the delay of 
    manually putting in the detection strip.
    This model includes the transformation of the graphs to allign the plateaus.

    Parameters:
    ----------
    parameters: np.ndarray
        Array containing the model parameters (float): transcription rate, translation rate, 
        maturation rate, catalytic rate of the enzyme, Michaelis constant transcription, 
        scaling factor transcription resources, Michaelis constant translation, Michaelis 
        constant translation resources, Michaelis constant enzymatic reaction, degradation 
        rate mRNA, degredation rate transcription resources association rate aptamer and 
        vitamin, dissociation rate aptamer and vitamin, cleaving rate aptamer.
    constants: np.ndarray
        Array containing the model constants (float): height of the microfluidic paper, 
        extinction coefficient CPRG, extinction coefficient CPR, blank intensity measurement 
        CPRG, blank intensity measurement CPR.
    initial_conditions: np.ndarray
        Array containing the initial conditions: gene concentration, substrate 
        concentration, vitamins.
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)

    Returns
    ----------
    time: np.ndarray
        Array containing all timepoints of the simulation.
    b_y: np.ndarray 
        Array containing the blue over yellow light intensity ratio at each timepoint.
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
    # I_null = constants[3] Has to be added after obtaining I0.

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
    # Concentration of beta-galactosidase (enzyme)
    e = np.zeros(n, dtype=np.float64)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    # Blue over yellow absorbance ratio
    b_y = np.zeros(n, dtype=np.float64)
    # Intentisity over time of the product Has to be added after obtaining I0. For now, A is used
    # I = np.zeros(n, dtype=np.float64)
    A = np.zeros(n, dtype=np.float64)

    # "Unpacking" the array with initial conditions into individual initial conditions
    dna[0] = initial_conditions[0]
    s[0] = initial_conditions[1]
    tsr[0] = 1
    tlr[0] = 1
    vit[0] = initial_conditions[2]

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
            0.25 * k_mat * e_mon[step]
        e_dt = 0.25 * k_mat * e_mon[step]
        s_dt = - k_cat * e[step] * s[step] / (k_m + s[step])
        p_dt = - s_dt
        # I_dt = I_null / ( 10**(p_dt * eps_cpr * h ) ) Has to be added after obtaining I0. For now, A is used.
        A_dt = p_dt*eps_cpr*h

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
        # I[step + 1] = I[step] + I_dt * dt Has to be added after obtaining I0. For now, A is used.
        A[step + 1] = A[step] + A_dt * dt

    # Calculating blue over yellow ratio
    blue = eps_cpr * p * h
    yellow = eps_cprg * s * h
    b_y = np.divide(blue, yellow)
    p = p / np.amax(p)

    return time, A


# if __name__ == "__main__":
#     parameters = standard_parameters_prokaryotic()
#     constants = standard_constants()
#     initial_conditions = standard_initial_conditions(
#         dna_conc=5*10**-3, s_i=150, vit_conc=5)

#     time, p = model_prokaryotic_hp(
#         parameters, constants, initial_conditions, dt=1, t_tot=7200, dt_store=0.01)
