import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from models import model_prokaryotic
from standard_parameters import standard_parameters


def micromolar_conc_to_math_exp(conc: float):
    """Function sets a float value that with micromolar units into a correct scientific math expression"""
    count = 0
    temp_conc = conc  # temporary variable to store a possible
    while temp_conc < 1:
        temp_conc = temp_conc * 10
        count += 1
    if count == 0:
        math_exp = f"${temp_conc:.2f}$ $\mu M$"
    else:
        math_exp = f"${temp_conc:.2f} \\times 10^{{-{count}}}$ $\mu M$"
    return math_exp


def anim_two_vitamin_conc_differing_dna_conc(vit_conc1, vit_conc2, low_dna_conc=1*10**-6, high_dna_conc=5*10**-3, num_steps=10, dt=0.01, t_tot=7200, save=False):
    """All inputs are in micromolar"""
    # Determine all the DNA concentrations to try
    dna_conc_all = np.linspace(low_dna_conc, high_dna_conc, num_steps)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    absorbance1 = np.zeros((num_steps, timesteps), dtype=np.float32)
    absorbance2 = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_step = dna_conc_all
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    for i in range(num_steps):
        parameters, constants, initial_conditions = standard_parameters(
            dna_conc=dna_conc_all[i], vit_conc=vit_conc1)
        time1, absorbance1[i, :] = model_prokaryotic(
            parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

        parameters, constants, initial_conditions = standard_parameters(
            dna_conc=dna_conc_all[i], vit_conc=vit_conc2)
        time2, absorbance2[i, :] = model_prokaryotic(
            parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

        area = np.sum((absorbance1[i, :] - absorbance2[i, :]) * dt)

        # area = np.sum(
        #     (np.cumsum(absorbance1[i, :]) - np.cumsum(absorbance2[i, :])) * dt)
        area_array[i] = area

    # absorbance1 = np.cumsum(absorbance1, axis=1)
    # absorbance2 = np.cumsum(absorbance2, axis=1)

    area_array = area_array / area_array[0]

    # Create the figure
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()
    ax2.plot(area_step, area_array)

    # Store the label str expressions
    label_line1 = micromolar_conc_to_math_exp(vit_conc1) + " original"
    label_line2 = micromolar_conc_to_math_exp(vit_conc2) + " original"
    label_line3 = micromolar_conc_to_math_exp(vit_conc1) + " moving"
    label_line4 = micromolar_conc_to_math_exp(vit_conc2) + " moving"

    line1, = ax.plot(time1, absorbance1[0, :],
                     label=label_line1, color="#E389BB")
    line3, = ax.plot(time1, absorbance1[0, :],
                     label=label_line3, color="#9B0138")
    line2, = ax.plot(time2, absorbance2[0, :],
                     label=label_line2, color="#8B992F")
    line4, = ax.plot(time2, absorbance2[0, :],
                     label=label_line4, color="#667817")

    dna_conc_math_exp = "DNA concentration:\n" + \
        micromolar_conc_to_math_exp(high_dna_conc)
    dna_conc_text = ax.text(100, 100, dna_conc_math_exp,
                            fontsize=10, bbox=dict(facecolor="#FFCF39", alpha=0.5, boxstyle="round"))

    def init():
        # ax.grid(True)
        ax.legend()
        ax.set_xlim(0, t_tot)
        # ax.set_ylim(-10, 160)
        return line1, line2, line3, line4, dna_conc_text,

    def update(index):
        dna_conc_math_exp = "DNA concentration:\n" + \
            micromolar_conc_to_math_exp(dna_conc_all[index])
        dna_conc_text.set_text(dna_conc_math_exp)

        line3.set_data(time1, absorbance1[index, :])
        line4.set_data(time2, absorbance2[index, :])
        return line1, line2, line3, line4, dna_conc_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=True)

    if save:
        f = "test.gif"
        writergif = PillowWriter(fps=60)
        anim.save(f, writer=writergif)

    plt.show()


anim_two_vitamin_conc_differing_dna_conc(
    0.5, 2, low_dna_conc=0.5*10**-4, high_dna_conc=5*10**-3, num_steps=50, dt=0.01, t_tot=14400)


# parameters, constants, initial_conditions = standard_parameters(
#     dna_conc=5*10**-3, vit_conc=5)
# time1, absorbance1 = model_prokaryotic(
#     parameters, constants, initial_conditions, dt=0.01, t_tot=7200)


# parameters, constants, initial_conditions = standard_parameters(
#     dna_conc=5*10**-3, vit_conc=20)
# time2, absorbance2 = model_prokaryotic(
#     parameters, constants, initial_conditions, dt=0.01, t_tot=7200)

# time_text = ax.text(100, 100,
#                     r"DNA concentration: $5 \times 10^{-3}$ $\mu M$", fontsize=10)

# graph_step = 50
# label_line1 = micromolar_conc_to_math_exp(5) + " original"
# line1, = plt.plot(time1[::graph_step],
#                   absorbance1[::graph_step], label=label_line1)
# line3, = plt.plot(time1[::graph_step],
#                   absorbance1[::graph_step], label="5 move")
# line2, = plt.plot(time2[::graph_step],
#                   absorbance2[::graph_step], label="50 standard")
# line4, = plt.plot(time2[::graph_step],
#                   absorbance2[::graph_step], label="50 move")

# n = 20
# timesteps = int(np.ceil(T_TOT/DT)) + 1
# absorbance1 = np.zeros((n, timesteps), dtype=np.float32)
# absorbance2 = np.zeros((n, timesteps), dtype=np.float32)


# dna_concenctrations = np.linspace(1*10**-6, 5*10**-3, n)[::-1]

# for i in range(n):
#     parameters, constants, initial_conditions = standard_parameters(
#         dna_conc=dna_concenctrations[i], vit_conc=5)
#     _, absorbance1[i, :] = model_prokaryotic(
#         parameters, constants, initial_conditions, dt=DT, t_tot=T_TOT)

#     parameters, constants, initial_conditions = standard_parameters(
#         dna_conc=dna_concenctrations[i], vit_conc=20)
#     _, absorbance2[i, :] = model_prokaryotic(
#         parameters, constants, initial_conditions, dt=DT, t_tot=T_TOT)


# def init():
#     ax.legend()
#     ax.set_xlim(0, 7200)
#     # ax.set_ylim(0.0, 2.0)
#     return line1, line2, line3, line4, time_text,


# def update(index):
#     label_dna = "DNA concentration: " + \
#         micromolar_conc_to_math_exp(dna_concenctrations[index])
#     time_text.set_text(label_dna)
#     line3.set_data(time1[::graph_step], absorbance1[index, ::graph_step])
#     line4.set_data(time2[::graph_step], absorbance2[index, ::graph_step])
#     return line1, line2, line3, line4, time_text,


# ani = FuncAnimation(fig, update, frames=np.arange(n),
#                     init_func=init, blit=True)

# plt.show()
