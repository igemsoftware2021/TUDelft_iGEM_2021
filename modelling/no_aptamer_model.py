# @njit
def no_aptamer_model(parameters, dt, t_tot, dna_i, s_i):
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
    # Concentration of mRNA [μM]
    mrna = np.zeros(n, dtype=np.float64)
    # Concentration of degraged mRNA [μM]
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

    # Plug in the initial concentrations
    dna[0] = dna_i
    s[0] = s_i
    tsr[0] = 1
    tlr[0] = 1

    # A loop with the differential equations
    for step in range(len(time)-1):
        # Differential of each species w.r.t time
        dna_dt = 0  # could remove this one, zero anyway
        mrna_dt = k_ts * tsr[step] * dna[step] / \
            (k_s + dna[step]) - deg_mrna * mrna[step]
        dmrna_dt = deg_mrna * mrna[step]
        tsr_dt = - kc_s * tsr[step] * dna[step] / (k_s + dna[step])
        tlr_dt = - deg_tlr * tlr[step] / (k_tlr + tlr[step])
        e_mon_dt = k_tl * tlr[step] * \
            mrna[step] / (k_l + mrna[step]) - \
            0.25 * k_mat * e_mon[step]
        e_dt = 0.25 * k_mat * e_mon[step]
        s_dt = - k_cat * e[step] * s[step] / (k_m + s[step])
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        dna[step + 1] = dna[step] + dna_dt * dt
        mrna[step + 1] = mrna[step] + mrna_dt * dt
        dmrna[step + 1] = dmrna[step] + dmrna_dt * dt
        tsr[step + 1] = tsr[step] + tsr_dt * dt
        tlr[step + 1] = tlr[step] + tlr_dt * dt
        e_mon[step + 1] = e_mon[step] + e_mon_dt * dt
        e[step + 1] = e[step] + e_dt * dt
        s[step + 1] = s[step] + s_dt * dt
        p[step + 1] = p[step] + p_dt * dt
    return (p, time)
