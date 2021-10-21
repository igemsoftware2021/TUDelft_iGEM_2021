# Libraries to import
import numpy as np
import matplotlib.pyplot as plt
from models import model_prokaryotic_absorbance
from standard_values import standard_parameters_prokaryotic, standard_constants
import alignment_helpers
from simulation_helpers import create_absorbance_sd, simulate_hardware

# Parameters
parameters = standard_parameters_prokaryotic()
# Constants
constants = standard_constants()

t_tot = 7200  # total time [s]
dt = 0.01  # timestep [s]

num_graph = 10  # Amount of samples to be calibrated with

# Empty array to store A(t) values for later
n = int(np.ceil(t_tot/dt)) + 1  # Number of timesteps of the simulation [-]
# Absorbance = []

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

absorbance_values, sd_values = create_absorbance_sd(20)

for i in np.linspace(5*10**-3, 5000*10**-3, num_graph):

    # Initial concentration of the beta-galactosidase gene [μM]
    dna_i = 3*10**-3
    s_i = 1000         # Initial substrate concentration [μM]

    initial_conditions = np.array([dna_i, s_i, i])

    # Running the model
    (time, absorbance) = model_prokaryotic_absorbance(parameters,
                                                      constants, initial_conditions, dt=dt, t_tot=t_tot)

    # time_hardware, absorbance_hardware = simulate_hardware(
    #     time, absorbance, absorbance_values, sd_values)

    absorbance = absorbance[::sample_size]
    time = time[::sample_size]
    time_array[:, count] = time
    absorbance_array[:, count] = absorbance
    count += 1

    ax1.plot(time, absorbance)
    # ax1.scatter(time_hardware, absorbance_hardware, s=0.5)

ax1.set_ylabel("Absorbance A [-]")
ax1.set_xlabel("Time t [s]")

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

plt.figure(12)
plt.imshow(a_rmsd_i)
plt.show()
