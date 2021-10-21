# Libraries to import
from models import model_prokaryotic
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from standard_values import standard_parameters_prokaryotic

parameters = standard_parameters_prokaryotic()

t_tot = 7200  # total time [s]
dt = 0.01  # timestep [s]

# Initial concentrations
dna_i = 3*10**-3  # Initial concentration of the beta-galactosidase gene [μM]
s_i = 1000  # Initial substrate concentration [μM]
vit_i = 0.07  # 0**-3  # Initial vitamin concentration [μM]
# Array containing above constants
initial_conditions = np.array([dna_i, s_i, vit_i])

# Running the model
(time, data) = model_prokaryotic(
    parameters, initial_conditions, dt=dt, t_tot=t_tot)

# Create parameter dictionairy
parameter_names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
                           "k_tlr", "k_m", "deg_mrna", "deg_tlr", "k_on", "k_off", "k_c", "vit_i", "dna_i", "s_i"]
parameters = np.append(parameters, vit_i)
parameters = np.append(parameters, dna_i)
parameters = np.append(parameters, s_i)

# Create the figure and the line
fig, ax = plt.subplots(figsize=(12, 8), dpi=150)
line, = plt.plot(time, data, lw=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel("Product concentration µM")
# ax.set_ylim(0, 1100)
axcolor = 'lightgoldenrodyellow'
ax.margins(x=0)

# adjust the main plot to make room for the sliders
# probably need more room
left = 0.25
plt.subplots_adjust(left=left, bottom=0.65)

# Make an array with axes that are sufficiently spaced
sliders_position = np.zeros([len(parameter_names), 4])
spacer = 0.03
for ii in range(len(parameter_names)):
    sliders_position[ii, 0] = left  # left
    sliders_position[ii, 1] = 0.55 - ii * spacer
    sliders_position[ii, 2] = 0.65  # width
    sliders_position[ii, 3] = 0.03  # height


axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
# Make a horizontal sliders to control each parameter.
parameter_sliders = []
for ii in range(len(parameter_names)):
    parameter_sliders.append(Slider(
        ax=plt.axes(sliders_position[ii, :], facecolor=axcolor),
        label=parameter_names[ii],
        valmin=(1/3) * parameters[ii],
        valmax=3 * parameters[ii],
        valinit=parameters[ii],
    ))

# The functions to be called anytime a slider's value changes


def update_k_ts(val):
    parameters[0] = parameter_sliders[0].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_tl(val):
    parameters[1] = parameter_sliders[1].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_mat(val):
    parameters[2] = parameter_sliders[2].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_cat(val):
    parameters[3] = parameter_sliders[3].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_s(val):
    parameters[4] = parameter_sliders[4].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_kc_s(val):
    slide_parameter = "kc_s"
    parameters[5] = parameter_sliders[5].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_l(val):
    parameters[6] = parameter_sliders[6].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_tlr(val):
    parameters[7] = parameter_sliders[7].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_m(val):
    parameters[8] = parameter_sliders[8].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_deg_mrna(val):
    parameters[9] = parameter_sliders[9].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_deg_tlr(val):
    parameters[10] = parameter_sliders[10].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_on(val):
    parameters[11] = parameter_sliders[11].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_off(val):
    parameters[12] = parameter_sliders[12].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_k_c(val):
    parameters[13] = parameter_sliders[13].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_vit_i(val):
    initial_conditions[2] = parameter_sliders[14].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_dna_i(val):
    initial_conditions[0] = parameter_sliders[15].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
    line.set_ydata(data)
    fig.canvas.draw_idle()


def update_s_i(val):
    initial_conditions[1] = parameter_sliders[16].val
    time, data = model_prokaryotic(parameters,
                                   initial_conditions, dt=dt, t_tot=t_tot)
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
parameter_sliders[11].on_changed(update_k_on)
parameter_sliders[12].on_changed(update_k_off)
parameter_sliders[13].on_changed(update_k_c)
parameter_sliders[14].on_changed(update_vit_i)
parameter_sliders[15].on_changed(update_dna_i)
parameter_sliders[16].on_changed(update_s_i)

plt.show()
