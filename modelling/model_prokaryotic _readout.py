# Libraries to import
from models import model_prokaryotic_readout
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import alignment_helpers

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
deg_tlr = 7.5*10**-5  # (10) Degradation rate of translation resources [1/s]
k_on = 1*10**-2       # (11)  Association rate of vitamin and umRNA [1/µMs]
k_off = 1*10**-2      # (12)  Dissociation rate of vitamin and umRNA [1/s]
k_c = (1/60)/10            # (13)  Cleaving rate of umRNA [1/s]
parameters = np.array([k_ts, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                       k_tlr, k_m, deg_mrna, deg_tlr, k_on, k_off, k_c])  # Array containing above parameters

t_tot = 7200  # total time [s]
dt = 0.01  # timestep [s]

num_graph = 10  # Amount of samples to be calibrated with

# Empty array to store A(t) values for later
n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
# Absorbance = []


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
sample_size = 20
time_array = np.zeros((int(np.ceil(n/sample_size)),
                      num_graph), dtype=np.float64)
absorbance_array = np.zeros(
    (int(np.ceil(n/sample_size)), num_graph), dtype=np.float64)
# Create the figure and the line
fig1, ax1 = plt.subplots()
count = 0
# For loop to create reference curves for now.
# vit concentration [µM]
for i in np.linspace(5*10**-3, 5000*10**-3, num_graph):

    # Initial concentration of the beta-galactosidase gene [μM]
    dna_i = 5*10**-3
    s_i = 150         # Initial substrate concentration [μM]

    initial_conditions = np.array([dna_i, s_i, i])

    # Running the model
    (time, absorbance) = model_prokaryotic_readout(parameters,
                                                   constants, initial_conditions, dt=dt, t_tot=t_tot)
    absorbance = absorbance[::sample_size]
    time = time[::sample_size]
    time_array[:, count] = time
    absorbance_array[:, count] = absorbance
    count += 1

    ax1.plot(time, absorbance)

ax1.set_ylabel("Absorbance A [-]")
ax1.set_xlabel("Time t [s]")
print(time_array.shape)

delta_x = np.array([0.01])
a_rsmd = np.zeros(len(delta_x), dtype=np.float64)
for i in range(len(delta_x)):
    x_shifted, y_shifted, a_rmsd_i = alignment_helpers.graph_alignment(
        time_array, absorbance_array, delta_x=delta_x[i], tol=1e-8)
    a_rsmd[i] = a_rmsd_i[-1][0]
a_rsmd = abs(a_rsmd-a_rsmd[-0])
# plt.figure(10)
#plt.loglog(delta_x, a_rsmd)

fig2, ax2 = plt.subplots()
for i in range(len(y_shifted)):
    ax2.plot(x_shifted, y_shifted[i])

fig1.show()
fig2.show()
plt.show()

# data = np.transpose(Absorbance)
# data = np.insert(data, 0, time, axis=1)

# df_data = pd.DataFrame(data)
# filepath = "calibration.xlsx"
# df_data.to_excel(filepath, index=False)

plt.figure(12)
plt.imshow(a_rmsd_i)
plt.show()
