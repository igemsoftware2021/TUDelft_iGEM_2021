import timeit
import numpy as np
from numba import njit
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions


def model_prokaryotic_python(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
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

    return time, p


@njit(cache=True, nogil=True)
def model_prokaryotic_numba(parameters, constants, initial_conditions, dt=0.01, t_tot=7200):
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

    return time, p


if __name__ == "__main__":
    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()
    initial_conditions = standard_initial_conditions(
        dna_conc=5*10**-3, s_i=150, vit_conc=1)

    t_python = timeit.Timer(
        "model_prokaryotic_python(parameters, constants, initial_conditions)", globals=globals())
    t_numba = timeit.Timer(
        "model_prokaryotic_numba(parameters, constants, initial_conditions)", globals=globals())

    number_python = 3
    number_numba = 3

    timings_python = t_python.repeat(repeat=5, number=number_python)
    timings_numba = t_numba.repeat(repeat=20, number=number_numba)

    min_timings_python = min(timings_python)
    min_timings_numba = min(timings_numba)

    min_timings_python_per_loop = min_timings_python/number_python
    min_timings_numba_per_loop = min_timings_numba/number_numba

    print(
        f"fastest loop python: {min_timings_python_per_loop}, fastest loop numba: {min_timings_numba_per_loop}")
    print(
        f"This is a speed increase of {int(min_timings_python_per_loop/min_timings_numba_per_loop)}x")
