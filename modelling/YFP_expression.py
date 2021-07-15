from models import model_yfp_expression
import numpy as np
import matplotlib.pyplot as plt
import csv

# Parameters
# (#) denotes the position in the parameters array
k_ts = 18.2 / 60 / 1000     # (0)  Transcription rate [uM/s]
k_tl = 7.2*10**-5     # (1)  Enzyme translation rate [uM/s]
k_mat = 0.5*10**-3    # (2)  Maturation rate of beta-galactosidase [1/s]
k_s = 8.5*10**-3      # (3)  Michaelis constant of transcription [μM]
# (4)  Scaling factor for the transcription resources [-]
kc_s = 1.8*10**-4
k_l = 65.8*10**-3     # (5)  Michaelis constant of translation [μM]
# (6) Michaelis constant of translation resources [-]
k_tlr = 6*10**-5
deg_mrna = 0  # (7) Degradation rate of mRNA [1/s]
# (8) Degradation rate of translation resources [1/s]
deg_tlr = 7.5*10**-5
parameters = np.array([k_ts, k_tl, k_mat, k_s, kc_s, k_l,
                       k_tlr, deg_mrna, deg_tlr])  # Array containing above parameters

t_tot = 60000  # total time [s]
dt = 0.01  # timestep [s]

# Initial DNA concentration
dna_i = 5*10**-3  # Initial concentration of the beta-galactosidase gene [μM]


# Running the model
time_model, yfp_model = model_yfp_expression(
    parameters, dna_i, dt=dt, t_tot=t_tot)


def cfs_reader(filename):
    file = open(filename, 'r')
    csv_reader = csv.reader(file)
    time = []
    data = []
    a = 1
    for line in csv_reader:
        if a == 1:
            time = [0]
            data = [261]
        else:
            time.append(np.float32(line[0].replace(",", ".")))
            data.append(np.float32(line[1].replace(",", ".")))
        a = a + 1
    time = np.array(time)
    data = np.array(data)
    data = data - data[0]
    return time, data


path = "modelling\data\cfs\\"
filename = "06_22_4YFP_Liquid_PURE2.csv"
time_cfs, yfp_cfs = cfs_reader(path + filename)


def normalize(data):
    data_max = np.amax(data)
    data_normalized = data / data_max
    return data_normalized


yfp_model = normalize(yfp_model)
yfp_cfs = normalize(yfp_cfs)

fig, ax = plt.subplots()
ax.plot(time_model, yfp_model, label="model", color="#9B0138")
ax.plot(time_cfs, yfp_cfs, label="cfs", color="#FFCF39")
ax.legend()
ax.set_xlabel("Time (s)"), ax.set_ylabel("YFP")
plt.show()
