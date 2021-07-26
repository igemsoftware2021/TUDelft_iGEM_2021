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

t_tot = 3 * 3600  # total time [s]
dt = 0.01  # timestep [s]

# Initial concentrations
dna_i = 5*10**-3  # Initial concentration of the beta-galactosidase gene [μM]
s_i = 1000   # Initial substrate concentration [μM]
# Array containing above constants
initial_conditions = np.array([dna_i, s_i])

# Constants
# (#) denotes the position in the constants array
h = 0.020  # (0) Height of the paper [cm]
eps_cprg = 0.294  # (1) Exctinction coefficient of CPRG at a wavelength of ???
eps_cpr = 0.539   # (2) Exctinction coefficient of CPR at a wavelength of ???
# Array containing above constants
constants = np.array([h, eps_cprg, eps_cpr])

# Running the model
time, data = model_no_aptamer(parameters,
                              constants, initial_conditions, dt=dt, t_tot=t_tot)


slide_parameter = "k_cat"


# Create parameter dictionairy
parameter_names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
                           "k_tlr", "k_m", "deg_mrna", "deg_tlr"]


# Create the figure and the line
fig, ax = plt.subplots()
line, = plt.plot(time, data, lw=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel("mRNA concentration uM")

axcolor = 'lightgoldenrodyellow'
ax.margins(x=0)

# adjust the main plot to make room for the sliders
# probably need more room
left = 0.25
plt.subplots_adjust(left=left, bottom=0.65)

# Make an array with axes that are sufficiently spaced
sliders_position = np.zeros([parameters.shape[0], 4])
spacer = 0.05
for ii in range(len(parameters)):
    sliders_position[ii, 0] = left  # left
    sliders_position[ii, 1] = 0.55 - ii * spacer
    sliders_position[ii, 2] = 0.65  # width
    sliders_position[ii, 3] = 0.03  # height


axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
# Make a horizontal sliders to control each parameter.
parameter_sliders = []
for ii in range(len(parameters)):
    parameter_sliders.append(Slider(
        ax=plt.axes(sliders_position[ii, :], facecolor=axcolor),
        label=parameter_names[ii],
        valmin=0.1 * parameters[ii],
        valmax=10 * parameters[ii],
        valinit=parameters[ii],
    ))

# The functions to be called anytime a slider's value changes


def update_k_ts(val):
    parameters[0] = parameter_sliders[0].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_tl(val):
    parameters[1] = parameter_sliders[1].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_mat(val):
    parameters[2] = parameter_sliders[2].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_cat(val):
    parameters[3] = parameter_sliders[3].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_s(val):
    parameters[4] = parameter_sliders[4].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_kc_s(val):
    slide_parameter = "kc_s"
    parameters[5] = parameter_sliders[5].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_l(val):
    parameters[6] = parameter_sliders[6].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_tlr(val):
    parameters[7] = parameter_sliders[7].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_m(val):
    parameters[8] = parameter_sliders[8].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_deg_mrna(val):
    parameters[9] = parameter_sliders[9].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_deg_tlr(val):
    parameters[10] = parameter_sliders[10].val
    time, data = model_no_aptamer(parameters,
                                  constants, initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


# register the update function with each slider

parameter_sliders[0].on_changed(update_k_ts)
parameter_sliders[1].on_changed(update_k_tl)
parameter_sliders[2].on_changed(update_k_mat)
parameter_sliders[3].on_changed(update_k_cat)
parameter_sliders[4].on_changed(update_k_s)
parameter_sliders[5].on_changed(update_kc_s)
parameter_sliders[6].on_changed(update_k_l)
parameter_sliders[7].on_changed(update_k_tlr)
parameter_sliders[8].on_changed(update_k_m)
parameter_sliders[9].on_changed(update_deg_mrna)
parameter_sliders[10].on_changed(update_deg_tlr)


# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
# resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
# button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


# def reset(event):
#   parameter_slider.reset()


# button.on_clicked(reset)

plt.show()
