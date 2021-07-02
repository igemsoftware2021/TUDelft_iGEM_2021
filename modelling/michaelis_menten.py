
import numpy as np
import matplotlib.pyplot as plt


def Michaelis_Menten(k_m, k_cat, dt, t_tot, s_i, e_i):
    # Determine the timepoints of the simulation
    n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints
    # Arrays for the concentration of each species)
    # Concentration of CPRG (substrate)
    s = np.zeros(n, dtype=np.float64)
    # Concentration of CPR (product)
    p = np.zeros(n, dtype=np.float64)
    e = np.zeros(n, dtype=np.float64)
    # Plug in the initial concentrations
    s[0] = s_i
    e[0] = e_i
    # A loop with the differential equations
    for step in range(len(time)-1):
        # Differential of each species w.r.t time
        s_dt = - k_cat * e[step] * s[step] / (k_m + s[step])
        p_dt = - s_dt

        # Computing the concentration of each species for the next step
        s[step + 1] = s[step] + s_dt * dt
        p[step + 1] = p[step] + p_dt * dt
        e[step + 1] = e[step]
    return s, p, e, time


k_m = 0.05  # [mM]
k_cat = 4280  # ong 4000 1/s
dt = 0.001
t_tot = 3600
s_i = 1  # [mM]
e = 1*10**-6  # [mM]
(s, p, e, time) = Michaelis_Menten(k_m, k_cat, dt, t_tot, s_i, e)
# Plotting to check if sth happened at all
plt.plot(time, p)
plt.show()
