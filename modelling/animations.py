import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegFileWriter
from models import model_prokaryotic, model_prokaryotic_readout, model_prokaryotic_all
from standard_values import standard_constants, standard_initial_conditions, standard_parameters_prokaryotic
from plot_helpers import micromolar_conc_to_math_exp


def anim_two_vitamin_conc_differing_dna_conc(vit_conc1, vit_conc2, low_dna_conc=1*10**-6, high_dna_conc=5*10**-3, num_steps=10, dt=0.01, t_tot=7200, save_path=None):
    """All inputs are in micromolar"""

    # The difference in dt between every plotted point, this is done to speed up the animation
    # so that if the dt=0.001 matplotlib does not have to plot every single point, while the
    # simulation can still be run at dt timesteps
    plot_dt = 0.1
    # The index difference in the timestep array to make sure you plot every point with
    # timestep difference of plot_dt. However the simulation is still run at timestep dt

    if dt >= plot_dt:
        plot_di = 1
    else:
        plot_di = int(np.floor((plot_dt / dt)))

    # Determine all the DNA concentrations to try
    dna_conc_all = np.linspace(low_dna_conc, high_dna_conc, num_steps)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    absorbance1 = np.zeros((num_steps, timesteps), dtype=np.float32)
    absorbance2 = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()

    for i in range(num_steps):
        initial_conditions = standard_initial_conditions(
            dna_conc=dna_conc_all[i], s_i=150, vit_conc=vit_conc1)
        time1, absorbance1[i, :] = model_prokaryotic_readout(
            parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

        initial_conditions = standard_initial_conditions(
            dna_conc=dna_conc_all[i], s_i=150, vit_conc=vit_conc2)
        time2, absorbance2[i, :] = model_prokaryotic_readout(
            parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

        area = np.sum((absorbance1[i, :] - absorbance2[i, :]) * dt)

        area_array[i] = area

    area_array = area_array / area_array[0]

    # Create the figure
    fig, (ax, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(
        12, 8), gridspec_kw={"height_ratios": [2, 1]})

    # fig2, ax2 = plt.subplots()

    # Store the label str expressions
    label_line1 = micromolar_conc_to_math_exp(vit_conc1, 0) + " original"
    label_line2 = micromolar_conc_to_math_exp(vit_conc2, 0) + " original"
    label_line3 = micromolar_conc_to_math_exp(vit_conc1, 0) + " moving"
    label_line4 = micromolar_conc_to_math_exp(vit_conc2, 0) + " moving"

    # First normal product lines
    line1, = ax.plot(time1[::plot_di], absorbance1[0, ::plot_di],
                     label=label_line1, color="#E389BB")
    line2, = ax.plot(time2[::plot_di], absorbance2[0, ::plot_di],
                     label=label_line2, color="#8B992F")
    line3, = ax.plot(time1[::plot_di], absorbance1[0, ::plot_di],
                     label=label_line3, color="#9B0138")
    line4, = ax.plot(time2[::plot_di], absorbance2[0, ::plot_di],
                     label=label_line4, color="#667817")

    fill1 = ax.fill_between(
        time1[::plot_di], absorbance1[0, ::plot_di], absorbance2[0, ::plot_di], color="#FFCF39", alpha=0.25)
    fill2 = ax.fill_between(
        time1[::plot_di], absorbance1[0, ::plot_di], absorbance2[0, ::plot_di], color="#FFCF39", alpha=0.75)

    # Area line
    area_line1, = ax2.plot(dna_conc_all, area_array, color="#E389BB")
    area_line2, = ax2.plot(dna_conc_all[0], area_array[0], color="#9B0138")
    # ax2.plot(dna_conc_all, area_array)

    dna_conc_math_exp = "DNA concentration: " + \
        micromolar_conc_to_math_exp(high_dna_conc, 2)
    dna_conc_text = ax.text(0.7, 0.1, dna_conc_math_exp, transform=ax.transAxes,
                            fontsize=10, bbox=dict(facecolor="#FFCF39", alpha=0.5, boxstyle="round"))

    def init():
        # ax.grid(True)

        ax.set_title("Absorbance over time")
        ax.legend()
        ax.set_xlabel(r"Time $[s]$")
        ax.set_ylabel(r"Absorbance $[AU]$")
        ax.set_xlim(0, t_tot)
        # ax.set_ylim(-10, 160)

        ax2.set_title("Relative area between two graphs over time")
        ax2.set_xlabel(r"DNA concentration $[\mu M]$")
        ax2.set_ylabel("Relative area")

        # Determine how to set the x-axis for ax2
        x_lim_diff = np.amax(dna_conc_all)*0.05
        ax2.set_xlim(np.amin(dna_conc_all)-x_lim_diff,
                     np.amax(dna_conc_all)+x_lim_diff)

        # Determine how to set the y-axis for ax2
        y_lim_diff = np.amax(area_array)*0.05
        ax2.set_ylim(np.amin(area_array)-y_lim_diff,
                     np.amax(area_array)+y_lim_diff)

        # Set a tight_layout for the figure. This needs to be done
        # after the axis names and title have been set
        fig.tight_layout()
        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, dna_conc_text,

    def update(index):

        dna_conc_math_exp = "DNA concentration: " + \
            micromolar_conc_to_math_exp(dna_conc_all[index], 2)
        dna_conc_text.set_text(dna_conc_math_exp)

        line3.set_data(time1[::plot_di], absorbance1[index, ::plot_di])
        line4.set_data(time2[::plot_di], absorbance2[index, ::plot_di])

        area_line2.set_data(dna_conc_all[:index], area_array[:index])

        ax.collections.clear()
        fill1 = ax.fill_between(
            time1[::plot_di], absorbance1[0, ::plot_di], absorbance2[0, ::plot_di], color="#FFCF39", alpha=0.25)
        fill2 = ax.fill_between(
            time1[::plot_di], absorbance1[index, ::plot_di], absorbance2[index, ::plot_di], color="#FFCF39", alpha=0.75)

        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, dna_conc_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=False)

    if save_path is not None:
        writermp4 = FFMpegFileWriter(fps=5, bitrate=5000)
        anim.save(f"{save_path}", writer=writermp4)

    plt.show()


def anim_two_vitamin_conc_differing_k_c(vit_conc1, vit_conc2, dna_conc, low_k_c=1*10**-6, high_k_c=5*10**-3, num_steps=10, dt=0.01, t_tot=7200, save=False):
    """All inputs are in micromolar"""

    # The difference in dt between every plotted point, this is done to speed up the animation
    # so that if the dt=0.001 matplotlib does not have to plot every single point, while the
    # simulation can still be run at dt timesteps
    plot_dt = 0.1
    # The index difference in the timestep array to make sure you plot every point with
    # timestep difference of plot_dt. However the simulation is still run at timestep dt

    if dt >= plot_dt:
        plot_di = 1
    else:
        plot_di = int(np.floor((plot_dt / dt)))

    # Determine all the DNA concentrations to try
    k_c_all = np.linspace(low_k_c, high_k_c, num_steps)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    absorbance1 = np.zeros((num_steps, timesteps), dtype=np.float32)
    absorbance2 = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()

    initial_conditions1 = standard_initial_conditions(
        dna_conc=dna_conc, s_i=150, vit_conc=vit_conc1)
    initial_conditions2 = standard_initial_conditions(
        dna_conc=dna_conc, s_i=150, vit_conc=vit_conc2)

    for i in range(num_steps):

        # Change the k_c parameter
        parameters[13] = k_c_all[i]

        time1, absorbance1[i, :] = model_prokaryotic_readout(
            parameters, constants, initial_conditions1, dt=dt, t_tot=t_tot)

        time2, absorbance2[i, :] = model_prokaryotic_readout(
            parameters, constants, initial_conditions2, dt=dt, t_tot=t_tot)

        area = np.sum((absorbance1[i, :] - absorbance2[i, :]) * dt)

        area_array[i] = area

    area_array = area_array / area_array[0]

    # Create the figure
    fig, (ax, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(
        12, 8), gridspec_kw={"height_ratios": [2, 1]})

    # fig2, ax2 = plt.subplots()

    # Store the label str expressions
    label_line1 = micromolar_conc_to_math_exp(vit_conc1, 0) + " original"
    label_line2 = micromolar_conc_to_math_exp(vit_conc2, 0) + " original"
    label_line3 = micromolar_conc_to_math_exp(vit_conc1, 0) + " moving"
    label_line4 = micromolar_conc_to_math_exp(vit_conc2, 0) + " moving"

    # First normal product lines
    line1, = ax.plot(time1[::plot_di], absorbance1[0, ::plot_di],
                     label=label_line1, color="#E389BB")
    line2, = ax.plot(time2[::plot_di], absorbance2[0, ::plot_di],
                     label=label_line2, color="#8B992F")
    line3, = ax.plot(time1[::plot_di], absorbance1[0, ::plot_di],
                     label=label_line3, color="#9B0138")
    line4, = ax.plot(time2[::plot_di], absorbance2[0, ::plot_di],
                     label=label_line4, color="#667817")

    fill1 = ax.fill_between(
        time1[::plot_di], absorbance1[0, ::plot_di], absorbance2[0, ::plot_di], color="#FFCF39", alpha=0.25)
    fill2 = ax.fill_between(
        time1[::plot_di], absorbance1[0, ::plot_di], absorbance2[0, ::plot_di], color="#FFCF39", alpha=0.75)

    # Area line
    area_line1, = ax2.plot(k_c_all, area_array, color="#E389BB")
    area_line2, = ax2.plot(k_c_all[0], area_array[0], color="#9B0138")
    # ax2.plot(dna_conc_all, area_array)

    k_c_math_exp = "Cleaving rate: " + \
        micromolar_conc_to_math_exp(high_k_c, 2)
    k_c_text = ax.text(0.7, 0.1, k_c_math_exp, transform=ax.transAxes,
                       fontsize=10, bbox=dict(facecolor="#FFCF39", alpha=0.5, boxstyle="round"))

    def init():
        # ax.grid(True)

        ax.set_title("Absorbance over time")
        ax.legend()
        ax.set_xlabel(r"Time $[s]$")
        ax.set_ylabel(r"Absorbance $[AU]$")
        ax.set_xlim(0, t_tot)
        # ax.set_ylim(-10, 160)

        ax2.set_title("Relative area between two graphs over time")
        ax2.set_xlabel(r"Cleaving rate $[1/s]]$")
        ax2.set_ylabel("Relative area")

        # Determine how to set the x-axis for ax2
        x_lim_diff = np.amax(k_c_all)*0.05
        ax2.set_xlim(np.amin(k_c_all)-x_lim_diff,
                     np.amax(k_c_all)+x_lim_diff)

        # Determine how to set the y-axis for ax2
        y_lim_diff = np.amax(area_array)*0.05
        ax2.set_ylim(np.amin(area_array)-y_lim_diff,
                     np.amax(area_array)+y_lim_diff)

        # Set a tight_layout for the figure. This needs to be done
        # after the axis names and title have been set
        fig.tight_layout()
        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, k_c_text,

    def update(index):

        k_c_math_exp = "Cleaving rate: " + \
            micromolar_conc_to_math_exp(k_c_all[index], 2)
        k_c_text.set_text(k_c_math_exp)

        line3.set_data(time1[::plot_di], absorbance1[index, ::plot_di])
        line4.set_data(time2[::plot_di], absorbance2[index, ::plot_di])

        area_line2.set_data(k_c_all[:index], area_array[:index])

        ax.collections.clear()
        fill1 = ax.fill_between(
            time1[::plot_di], absorbance1[0, ::plot_di], absorbance2[0, ::plot_di], color="#FFCF39", alpha=0.25)
        fill2 = ax.fill_between(
            time1[::plot_di], absorbance1[index, ::plot_di], absorbance2[index, ::plot_di], color="#FFCF39", alpha=0.75)

        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, k_c_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=False)

    if save:
        f = "test.mp4"
        writermp4 = FFMpegFileWriter(fps=5, bitrate=5000)
        anim.save(f, writer=writermp4)

    plt.show()


def anim_frac_mrna_conc_differing_dna_conc(vit_conc1, low_dna_conc=1*10**-6, high_dna_conc=5*10**-3, num_steps=10, dt=0.01, t_tot=7200, save=False):
    # TODO fix this function, this still uses old functions and incorrect names
    """All inputs are in micromolar"""
    # Determine all the DNA concentrations to try
    dna_conc_all = np.linspace(low_dna_conc, high_dna_conc, num_steps)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    frac_umrna = np.zeros((num_steps, timesteps), dtype=np.float32)
    frac_umrna_vit = np.zeros((num_steps, timesteps), dtype=np.float32)
    frac_cmrna = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_step = dna_conc_all
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    for i in range(num_steps):
        parameters, constants, initial_conditions = standard_parameters(
            dna_conc=dna_conc_all[i], vit_conc=vit_conc1)
        results = model_prokaryotic_all(
            parameters, constants, initial_conditions, dt=dt, t_tot=t_tot)

        if i == 0:
            time1 = results[0]

        umrna = results[2]
        umrna_vit = results[6]
        cmrna = results[3]
        total = umrna + umrna_vit + cmrna

        frac_umrna[i, :] = umrna / total
        frac_umrna_vit[i, :] = umrna_vit / total
        frac_cmrna[i, :] = cmrna / total

        # area = np.sum((absorbance1[i, :] - absorbance2[i, :]) * dt)

        # area = np.sum(
        #     (np.cumsum(absorbance1[i, :]) - np.cumsum(absorbance2[i, :])) * dt)
        # area_array[i] = area

    # absorbance1 = np.cumsum(absorbance1, axis=1)
    # absorbance2 = np.cumsum(absorbance2, axis=1)

    # area_array = area_array / area_array[0]

    # Create the figure
    fig, ax = plt.subplots()
    # fig2, ax2 = plt.subplots()
    # ax2.plot(area_step, area_array)

    # Store the label str expressions
    label_line1 = "umRNA"
    label_line2 = "umRNA_vit"
    label_line3 = "cmrna"

    line1, = ax.plot(time1, frac_umrna[0, :],
                     label=label_line1, color="#E389BB")
    line2, = ax.plot(time1, frac_umrna_vit[0, :],
                     label=label_line2, color="#8B992F")
    line3, = ax.plot(time1, frac_cmrna[0, :],
                     label=label_line3, color="#9B0138")
    # line4, = ax.plot(time2, absorbance2[0, :],
    #                  label=label_line4, color="#667817")

    dna_conc_math_exp = "DNA concentration:\n" + \
        micromolar_conc_to_math_exp(high_dna_conc, 2)
    dna_conc_text = ax.text(600, 0.8, dna_conc_math_exp,
                            fontsize=10, bbox=dict(facecolor="#FFCF39", alpha=0.5, boxstyle="round"))

    def init():
        # ax.grid(True)
        ax.legend()
        ax.set_xlim(0, t_tot)
        # ax.set_ylim(-10, 160)
        return line1, line2, line3, dna_conc_text,

    def update(index):
        dna_conc_math_exp = "DNA concentration:\n" + \
            micromolar_conc_to_math_exp(dna_conc_all[index], 2)
        dna_conc_text.set_text(dna_conc_math_exp)

        line1.set_data(time1, frac_umrna[index, :])
        line2.set_data(time1, frac_umrna_vit[index, :])
        line3.set_data(time1, frac_cmrna[index, :])

        return line1, line2, line3, dna_conc_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=True)

    if save:
        f = "test.gif"
        writergif = FFMpegFileWriter(fps=60)
        anim.save(f, writer=writergif)

    plt.show()


if __name__ == "__main__":
    # anim_two_vitamin_conc_differing_dna_conc(
    #     0.5, 1, low_dna_conc=0.5*10**-4, high_dna_conc=5*10**-3, num_steps=50, dt=0.01, t_tot=7200, save=False)

    # anim_two_vitamin_conc_differing_k_c(6, 7, 2*10**-3, low_k_c=(
    #     1/60)/10, high_k_c=(1/60), num_steps=30, dt=0.01, t_tot=7200, save=False)
