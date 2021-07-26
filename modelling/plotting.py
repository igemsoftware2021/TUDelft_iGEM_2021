import matplotlib.pyplot as plt
import numpy as np
from morris_method import morris_datareader
from models import model_no_aptamer

# Parameters
# (#) denotes the position in the parameters array
k_ts = 6.7*10**-5     # (0)  Transcription rate [uM/s]
k_tl = 7.2*10**-5     # (1)  Enzyme translation rate [uM/s]
k_mat = 0.5*10**-3    # (2)  Maturation rate of beta-galactosidase [1/s]
k_cat = 5.14*10**1   # (3)  Catalytic rate of beta-galactosidase [1/s]
k_s = 8.5*10**-3      # (4)  Michaelis constant of transcription [μM]
# (5)  Scaling factor for the transcription resources [-]
kc_s = 1.8*10**-4
k_l = 65.8*10**-3     # (6)  Michaelis constant of translation [μM]
# (7) Michaelis constant of translation resources [-]
k_tlr = 6*10**-5
k_m = 50               # (8) Michaelis constant of beta-galactosidase [μM]
deg_mrna = 1.3*10**-5  # (9) Degradation rate of mRNA [1/s]
# (10) Degradation rate of translation resources [1/s]
deg_tlr = 7.5*10**-5
parameters = np.array([k_ts, k_tl, k_mat, k_cat, k_s, kc_s, k_l,
                       k_tlr, k_m, deg_mrna, deg_tlr])  # Array containing above parameters

t_tot = 3600  # total time [s]
dt = 0.01  # timestep [s]

# Initial concentrations
dna_i = 5*10**-3  # Initial concentration of the beta-galactosidase gene [μM]
s_i = 1   # Initial substrate concentration
# Array containing above constants
initial_conditions = np.array([dna_i, s_i])

# Constants
# (#) denotes the position in the constants array
h = 0.02  # (0) Height of the paper [cm]
eps_cprg = 0.294  # (1) Exctinction coefficient of CPRG at a wavelength of ???
eps_cpr = 0.539   # (2) Exctinction coefficient of CPR at a wavelength of ???
# Array containing above constants
constants = np.array([h, eps_cprg, eps_cpr])

# Running the model
time, data = model_no_aptamer(parameters,
                              constants, initial_conditions, dt=dt, t_tot=t_tot)

# Plot every parameter of a single run

filenumber = "6"
tag = "_"
names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
         "k_tlr", "k_m", "deg_mrna", "deg_tlr"]
path = "modelling\data\morris_no_aptamer"


fig1, ax1 = plt.subplots()
time = morris_datareader("time", "mu", path, filenumber)
for name in names:
    mu_star = morris_datareader(name, "mu_star", path, filenumber)
    ax1.plot(time, mu_star, label="mu_star")

#ax1.set_title("Varied 1 order of magnitude")
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax1.legend(names, bbox_to_anchor=(1, 0.5), loc='center left',)
ax1.set_xlabel("Time (s)"), ax1.set_ylabel("Product (mM)")
plt.show()
fig1.savefig("modelling\data\morris_no_aptamer\\" +
             "all_6_mu_star" ".svg", format="svg", dpi=1200)


# Compare from 2 runs

# filenumber_1 = "5"
# filenumber_2 = "6"
# tag = "_"
# names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
#          "k_tlr", "k_m", "deg_mrna", "deg_tlr"]
# path = "modelling\data\morris_no_aptamer"


# time = morris_datareader("time", "mu", path, filenumber_1)
# for name in names:
#     fig1, ax1 = plt.subplots()
#     mu_star_1 = morris_datareader(name, "mu_star", path, filenumber_1)
#     mu_star_2 = morris_datareader(name, "mu_star", path, filenumber_2)
#     sigma_1 = morris_datareader(name, "sigma", path, filenumber_1)
#     sigma_2 = morris_datareader(name, "sigma", path, filenumber_2)
#     mu_1 = morris_datareader(name, "mu", path, filenumber_1)
#     mu_2 = morris_datareader(name, "mu", path, filenumber_2)
#     mu_star_conf_level_1 = morris_datareader(
#        name, "mu_star_conf_level", path, filenumber_1)
#     mu_star_conf_level_2 = morris_datareader(
#        name, "mu_star_conf_level", path, filenumber_2)
#     ax1.plot(time, mu_star_1/1000, label="mu_star 1 " + name, color="#E389BB")
#     ax1.plot(time, mu_star_2/1000, label="mu_star 3 " + name, color="#9B0138")
#     ax1.plot(time, sigma_1/1000, label="sigma 1 " + name, color="#8B992F")
#     ax1.plot(time, sigma_2/1000, label="sigma 3 " + name, color="#006400")
#     ax1.set_title(name)
#     ax1.set_ylim([0, 1])
#     ax1.legend()
#     ax1.set_xlabel("Time (s)"), ax1.set_ylabel("Product (mM)")
#     fig1.savefig("modelling\data\morris_no_aptamer\\" + name +
#                  "_combined_" + filenumber_1 + "_" + filenumber_2 + ".svg", format="svg", dpi=1200)
# plt.show()


# Plot every parameter of a single run


# filenumber = "3"
# tag = "_2"
# names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
#          "k_tlr", "k_m", "deg_mrna", "deg_tlr"]
# path = "modelling\data\morris_no_aptamer"


# for name in names:
#     fig1, ax1 = plt.subplots()
#     # ax1.plot(time, data, label="product", color="#9B0138")
#     mu_star = morris_datareader(name, "mu_star", path, filenumber)
#     #sigma = morris_datareader(name, "sigma", path, filenumber)
#     #mu = morris_datareader(name, "mu", path, filenumber)
#     # mu_star_conf_level = morris_datareader(
#     #    name, "mu_star_conf_level", path, filenumber)
#     # ax1.plot(time, mu, label="mu", color="#E389BB")
#     ax1.plot(time, mu_star, label="mu_star", color="#FFCF39")
#     # ax1.plot(time, sigma, label="sigma", color="#667817")
#     ax1.legend()
#     ax1.set_title(name)
#     ax1.set_xlabel("Time (s)"), ax1.set_ylabel("Product (mM)")
#     fig1.savefig("modelling\data\morris_no_aptamer\\" + name +
#                  "_" + filenumber + tag + ".svg", format="svg", dpi=1200)

# plt.show()
