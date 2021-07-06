# Libraries to import
from models import model_no_aptamer
from morris_method import morris_datareader
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import csv

tbd = 1  # Placeholder for unkown parameters/concentrations so that they're easy to find

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
h = 8*10**-5 * 1  # (0) Height of the paper [cm]
eps_cprg = 0.294  # (1) Exctinction coefficient of CPRG at a wavelength of ???
eps_cpr = 0.539   # (2) Exctinction coefficient of CPR at a wavelength of ???
# Array containing above constants
constants = np.array([h, eps_cprg, eps_cpr])

# Running the model
(time, data) = model_no_aptamer(parameters,
                                constants, initial_conditions, dt=dt, t_tot=t_tot)


slide_parameter = "k_ts"


# Create parameter dictionairy
parameter_dict = dict(zip(["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
                           "k_tlr", "k_m", "deg_mrna", "deg_tlr"], range(len(parameters))))
# Create the figure and the line
fig, ax = plt.subplots()
line, = plt.plot(time, data, lw=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel("Product concentration uM")

axcolor = 'lightgoldenrodyellow'
ax.margins(x=0)

# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.25, bottom=0.25)
# Make a horizontal slider to control the frequency.
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
parameter_slider = Slider(
    ax=axfreq,
    label=slide_parameter,
    valmin=0.01 * parameters[parameter_dict[slide_parameter]],
    valmax=10 * parameters[parameter_dict[slide_parameter]],
    valinit=parameters[parameter_dict[slide_parameter]],
)

# The function to be called anytime a slider's value changes


def update(val):
    parameters[parameter_dict[slide_parameter]] = parameter_slider.val
    (time, data) = model_no_aptamer(parameters,
                                    constants, initial_conditions, dt=dt, t_tot=7200)
    line.set_ydata(data)
    fig.canvas.draw_idle()


# register the update function with each slider
parameter_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    parameter_slider.reset()


button.on_clicked(reset)

plt.show()

##
# names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
#          "k_tlr", "k_m", "deg_mrna", "deg_tlr"]
# path = "modelling\data\morris_no_aptamer"
# filenumber = "1"
# time = morris_datareader("time", "mu", path, filenumber)
# for name in names:
#     mu_star = morris_datareader(name, "mu_star", path, filenumber)
#     sigma = morris_datareader(name, "sigma", path, filenumber)
#     mu = morris_datareader(name, "mu", path, filenumber)
#     mu_star_conf_level = morris_datareader(
#         name, "mu_star_conf_level", path, filenumber)
#     fig1 = plt.figure()
#     plt.fill_between(time, mu_star - 0.5 * mu_star_conf_level, mu_star + 0.5 *
#                      mu_star_conf_level, color="#8B992F")
#     plt.plot(time, data, label="product", color="#9B0138")
#     plt.plot(time, mu, label="mu", color="#FFCF39")
#     plt.plot(time, mu_star, label="mu_star", color="#667817")
#     plt.plot(time, sigma, label="sigma", color="#E389BB")
#     plt.title(name)
#     plt.legend()
#     plt.xlabel("Time (s)"), plt.ylabel("Product (uM)")
#     fig1.savefig("modelling\data\morris_no_aptamer\\" + name +
#                  "_" + filenumber + "_2" ".svg", format="svg", dpi=1200)

# plt.show()
