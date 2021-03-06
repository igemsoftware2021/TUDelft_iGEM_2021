import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegFileWriter
from models import model_prokaryotic, model_prokaryotic_all
from standard_values import standard_parameters_prokaryotic
from plot_helpers import micromolar_conc_to_math_exp
from matplotlib.ticker import MultipleLocator


def anim_two_vitamin_conc_differing_dna_conc(vit_conc1, vit_conc2, s_i=1000, low_dna_conc=1*10**-6, standard_dna_conc=3*10**-3, high_dna_conc=5*10**-3, num_steps=10, dt=0.01, t_tot=7200, save_path=None):
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
    dna_conc_all = np.geomspace(
        low_dna_conc, high_dna_conc, num=num_steps, endpoint=True)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    p1 = np.zeros((num_steps, timesteps), dtype=np.float32)
    p2 = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    parameters = standard_parameters_prokaryotic()

    # First calculate the standard area
    initial_conditions1 = np.array([standard_dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([standard_dna_conc, s_i, vit_conc2])

    time1, p1_standard = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)
    time2, p2_standard = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)
    area_standard = np.sum((p1_standard - p2_standard) * dt)

    for i in range(num_steps):
        initial_conditions1 = np.array([dna_conc_all[i], s_i, vit_conc1])
        initial_conditions2 = np.array([dna_conc_all[i], s_i, vit_conc2])

        time1, p1[i, :] = model_prokaryotic(
            parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)

        time2, p2[i, :] = model_prokaryotic(
            parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)

        area = np.sum((p1[i, :] - p2[i, :]) * dt)

        area_array[i] = area

    area_array = area_array / area_standard

    # Create the figure
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(
        14, 10), gridspec_kw={"height_ratios": [2, 1]}, dpi=125)

    # Store the label str expressions
    label_line1 = f"${int(vit_conc1*1000)}\;\mathrm{{nM}}$ - $\mathrm{{DNA}}$ fixed at {int(standard_dna_conc*1000)} $\mathrm{{nM}}$"
    label_line2 = f"${int(vit_conc2*1000)}\;\mathrm{{nM}}$ - $\mathrm{{DNA}}$ fixed at {int(standard_dna_conc*1000)} $\mathrm{{nM}}$"
    label_line3 = f"${int(vit_conc1*1000)}\;\mathrm{{nM}}$ - $\mathrm{{DNA}}$ variable"
    label_line4 = f"${int(vit_conc2*1000)}\;\mathrm{{nM}}$ - $\mathrm{{DNA}}$ variable"

    # First normal product lines
    line1, = ax1.plot(time1[::plot_di], p1_standard[::plot_di],
                      label=label_line1, color="#F3758A")
    line2, = ax1.plot(time2[::plot_di], p2_standard[::plot_di],
                      label=label_line2, color="#4FD590")
    line3, = ax1.plot(time1[::plot_di], p1[0, ::plot_di],
                      label=label_line3, color="#9B0138")
    line4, = ax1.plot(time2[::plot_di], p2[0, ::plot_di],
                      label=label_line4, color="#057D54")

    fill1 = ax1.fill_between(
        time1[::plot_di], p1_standard[::plot_di], p2_standard[::plot_di], color="#FFCE3A", alpha=0.25)
    fill2 = ax1.fill_between(
        time1[::plot_di], p1[0, ::plot_di], p2[0, ::plot_di], color="#FFCE3A", alpha=0.75)

    # Simulations have been ran, so make it nM
    dna_conc_all_nM = dna_conc_all * 1000
    standard_dna_conc_nM = standard_dna_conc * 1000

    # Area line
    area_line1, = ax2.plot(dna_conc_all_nM, area_array, color="#F3758A")
    area_line2, = ax2.plot(dna_conc_all_nM[0], area_array[0], color="#9B0138")
    # Plot the position of the dna_conc that determines the standard area
    ax2.scatter(standard_dna_conc_nM, 1.0, color="#057D54")

    dna_conc_math_exp = f"DNA concentration: {(standard_dna_conc_nM):.3f} $\mathrm{{nM}}$"
    dna_conc_text = ax1.text(0.01, 0.7, dna_conc_math_exp, transform=ax1.transAxes,
                             fontsize=10, bbox=dict(facecolor="#FFCE3A", alpha=0.5, boxstyle="round"))

    def init():
        # ax.set_title("Product concentration over time")
        ax1.legend()
        ax1.set_xlabel(r"Time $[\mathrm{{s}}]$")
        ax1.set_ylabel(r"Product concentration $[\mathrm{{\mu M}}]$")
        ax1.set_xlim(0, t_tot)

        # ax2.set_title("Area fold change between two graphs over time")
        ax2.set_xlabel(r"DNA concentration $[\mathrm{{n M}}]$")
        ax2.set_ylabel("Area fold change")

        # Determine how to set the x-axis for ax2
        x_lim_diff = np.amax(dna_conc_all_nM)*0.05
        ax2.set_xlim(np.amin(dna_conc_all_nM)-x_lim_diff,
                     np.amax(dna_conc_all_nM)+x_lim_diff)

        # Determine how to set the y-axis for ax2
        y_lim_diff = np.amax(area_array)*0.05
        ax2.set_ylim(np.amin(area_array)-y_lim_diff,
                     np.amax(area_array)+y_lim_diff)

        # Plot vertical and horizontal lines to the the position of the
        # dna_conc that determines the standard area
        ax2.vlines(standard_dna_conc_nM, 0, 1,
                   color="#057D54", linestyle="dashed")
        ax2.axhline(1, color="#057D54", linestyle="dashed")

        # Set minor and major tick locator
        ax2.xaxis.set_major_locator(MultipleLocator(0.5))
        ax2.xaxis.set_minor_locator(MultipleLocator(0.1))

        # Set the character labels
        ax1.text(-0.05, 1.05, "a", transform=ax1.transAxes,
                 size=16, weight="bold")
        ax2.text(-0.05, 1.05, "b", transform=ax2.transAxes,
                 size=16, weight="bold")

        # Set a tight_layout for the figure. This needs to be done
        # after the axis names and title have been set
        # fig.tight_layout()
        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, dna_conc_text,

    def update(index):

        dna_conc_math_exp = f"DNA concentration: {(dna_conc_all[index]*1000):.3f} $\mathrm{{nM}}$"
        dna_conc_text.set_text(dna_conc_math_exp)

        line3.set_data(time1[::plot_di], p1[index, ::plot_di])
        line4.set_data(time2[::plot_di], p2[index, ::plot_di])

        area_line2.set_data(dna_conc_all_nM[:index], area_array[:index])

        ax1.collections.clear()
        fill1 = ax1.fill_between(
            time1[::plot_di], p1_standard[::plot_di], p2_standard[::plot_di], color="#FFCE3A", alpha=0.25)
        fill2 = ax1.fill_between(
            time1[::plot_di], p1[index, ::plot_di], p2[index, ::plot_di], color="#FFCE3A", alpha=0.75)

        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, dna_conc_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=False)

    if save_path is not None:
        writermp4 = FFMpegFileWriter(fps=24, bitrate=24000)
        anim.save(f"{save_path}", writer=writermp4)
    else:
        plt.show()


def anim_two_vitamin_conc_differing_k_c(vit_conc1, vit_conc2, s_i=1000, dna_conc=3*10**-3, low_k_c=(1/60)/10, standard_k_c=1/60, high_k_c=(1/60)*10, num_steps=10, dt=0.01, t_tot=7200, save_path=None):
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

    # Determine all the cleaving rates to try
    k_c_all = np.geomspace(
        low_k_c, high_k_c, num=num_steps, endpoint=True)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    p1 = np.zeros((num_steps, timesteps), dtype=np.float32)
    p2 = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    parameters = standard_parameters_prokaryotic()
    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])

    # First calculate the standard area
    parameters[13] = standard_k_c
    time1, p1_standard = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)
    time2, p2_standard = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)
    area_standard = np.sum((p1_standard - p2_standard) * dt)

    for i in range(num_steps):
        parameters[13] = k_c_all[i]

        time1, p1[i, :] = model_prokaryotic(
            parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)

        time2, p2[i, :] = model_prokaryotic(
            parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)

        area = np.sum((p1[i, :] - p2[i, :]) * dt)

        area_array[i] = area

    area_array = area_array / area_standard

    # Create the figure
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(
        12, 8), gridspec_kw={"height_ratios": [2, 1]}, dpi=150)

    # Store the label str expressions
    label_line1 = f"${int(vit_conc1*1000)}\;\mathrm{{nM}}$ - $k_\mathrm{{c}}$ fixed at $0.017\;\mathrm{{s}}^{{-1}}$"
    label_line2 = f"${int(vit_conc2*1000)}\;\mathrm{{nM}}$ - $k_\mathrm{{c}}$ fixed at $0.017\;\mathrm{{s}}^{{-1}}$"
    label_line3 = f"${int(vit_conc1*1000)}\;\mathrm{{nM}}$ - $k_\mathrm{{c}}$ variable"
    label_line4 = f"${int(vit_conc2*1000)}\;\mathrm{{nM}}$ - $k_\mathrm{{c}}$ variable"

    # First normal product lines
    line1, = ax1.plot(time1[::plot_di], p1_standard[::plot_di],
                      label=label_line1, color="#F3758A")
    line2, = ax1.plot(time2[::plot_di], p2_standard[::plot_di],
                      label=label_line2, color="#4FD590")
    line3, = ax1.plot(time1[::plot_di], p1[0, ::plot_di],
                      label=label_line3, color="#9B0138")
    line4, = ax1.plot(time2[::plot_di], p2[0, ::plot_di],
                      label=label_line4, color="#057D54")

    fill1 = ax1.fill_between(
        time1[::plot_di], p1_standard[::plot_di], p2_standard[::plot_di], color="#FFCE3A", alpha=0.25)
    fill2 = ax1.fill_between(
        time1[::plot_di], p1[0, ::plot_di], p2[0, ::plot_di], color="#FFCE3A", alpha=0.75)

    # Area line
    area_line1, = ax2.plot(k_c_all, area_array, color="#F3758A")
    area_line2, = ax2.plot(k_c_all[0], area_array[0], color="#9B0138")
    # Plot the position of the k_c that is the standard area
    ax2.scatter(standard_k_c, 1.0, color="#057D54")

    k_c_math_exp = f"$k_{{\mathrm{{c}}}} = {k_c_all[0]:.4f}\;\mathrm{{s}}^{{-1}}$"
    k_c_text = ax1.text(0.01, 0.65, k_c_math_exp, transform=ax1.transAxes,
                        fontsize=10, bbox=dict(facecolor="#FFCE3A", alpha=0.5, boxstyle="round"))

    def init():
        # ax1.set_title("Product concentration over time")
        # ax1.legend(title="Vitamin concentration")
        ax1.legend()
        ax1.set_xlabel(r"Time $[\mathrm{{s}}]$")
        ax1.set_ylabel(r"Product concentration $[\mathrm{{\mu M}}]$")
        ax1.set_xlim(0, t_tot)

        # ax2.set_title("Area fold change between two graphs over time")
        ax2.set_xlabel(r"$k_\mathrm{c}\;[\mathrm{{s}}^{{-1}}]$")
        ax2.set_ylabel("Fold change")

        # Determine how to set the x-axis for ax2
        x_lim_diff = np.amax(k_c_all)*0.05
        ax2.set_xlim(np.amin(k_c_all)-x_lim_diff,
                     np.amax(k_c_all)+x_lim_diff)

        # Determine how to set the y-axis for ax2
        y_lim_diff = np.amax(area_array)*0.05
        ax2.set_ylim(np.amin(area_array)-y_lim_diff,
                     np.amax(area_array)+y_lim_diff)
        ax2.yaxis.set_major_locator(MultipleLocator(0.5))

        # Plot vertical and horizontal lines to the the position of the
        # k_c that determines the standard area
        ax2.vlines(standard_k_c, 0, 1, color="#057D54", linestyle="dashed")
        ax2.axhline(1, color="#057D54", linestyle="dashed")
        # ax2.hlines(1, np.amin(k_c_all)-x_lim_diff, standard_k_c,
        #            color="#057D54", linestyle="dashed")

        # Set the character labels
        ax1.text(-0.05, 1.05, "a", transform=ax1.transAxes,
                 size=16, weight="bold")
        ax2.text(-0.05, 1.05, "b", transform=ax2.transAxes,
                 size=16, weight="bold")

        # Set minor and major tick locator
        ax2.xaxis.set_major_locator(MultipleLocator(0.005))
        ax2.xaxis.set_minor_locator(MultipleLocator(0.001))

        # Set a tight_layout for the figure. This needs to be done
        # after the axis names and title have been set
        # fig.tight_layout()
        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, k_c_text,

    def update(index):

        k_c_math_exp = f"$k_{{\mathrm{{c}}}} = {k_c_all[index]:.4f}\;\mathrm{{s}}^{{-1}}$"
        k_c_text.set_text(k_c_math_exp)

        line3.set_data(time1[::plot_di], p1[index, ::plot_di])
        line4.set_data(time2[::plot_di], p2[index, ::plot_di])

        area_line2.set_data(k_c_all[:index], area_array[:index])

        ax1.collections.clear()
        fill1 = ax1.fill_between(
            time1[::plot_di], p1_standard[::plot_di], p2_standard[::plot_di], color="#FFCE3A", alpha=0.25)
        fill2 = ax1.fill_between(
            time1[::plot_di], p1[index, ::plot_di], p2[index, ::plot_di], color="#FFCE3A", alpha=0.75)

        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, k_c_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=False)

    if save_path is not None:
        writermp4 = FFMpegFileWriter(fps=24, bitrate=24000)
        anim.save(f"{save_path}", writer=writermp4)
    else:
        plt.show()


def anim_two_vitamin_conc_differing_k_D(vit_conc1, vit_conc2, s_i=1000, dna_conc=3*10**-3, low_k_D=0.005, standard_k_D=0.05, high_k_D=5, num_steps=10, dt=0.01, t_tot=7200, save_path=None):
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

    # Determine all the cleaving rates to try
    k_D_all = np.geomspace(
        low_k_D, high_k_D, num=num_steps, endpoint=True)[::-1]

    # Preallocate all the necessary storage
    timesteps = int(np.ceil(t_tot/dt)) + 1
    p1 = np.zeros((num_steps, timesteps), dtype=np.float32)
    p2 = np.zeros((num_steps, timesteps), dtype=np.float32)

    # Preallocate necessary storage
    area_array = np.zeros(num_steps, dtype=np.float32)

    # Precompute everything
    parameters = standard_parameters_prokaryotic()
    k_on_normal = parameters[11]
    k_off_normal = parameters[12]
    k_D_normal = k_off_normal/k_on_normal   # The k_D from the standard parameters

    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])

    # First calculate the standard area
    scaling_factor = k_D_normal/standard_k_D
    k_on_temp = k_on_normal * np.sqrt(scaling_factor)
    k_off_temp = k_off_normal / np.sqrt(scaling_factor)

    parameters[11] = k_on_temp
    parameters[12] = k_off_temp
    time1, p1_standard = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)
    time2, p2_standard = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)
    area_standard = np.sum((p1_standard - p2_standard) * dt)

    for i in range(num_steps):
        scaling_factor = k_D_normal / k_D_all[i]
        k_on_temp = k_on_normal * np.sqrt(scaling_factor)
        k_off_temp = k_off_normal / np.sqrt(scaling_factor)

        parameters[11] = k_on_temp
        parameters[12] = k_off_temp

        time1, p1[i, :] = model_prokaryotic(
            parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)

        time2, p2[i, :] = model_prokaryotic(
            parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)

        area = np.sum((p1[i, :] - p2[i, :]) * dt)

        area_array[i] = area

    area_array = area_array / area_standard

    # Create the figure
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(
        24, 10), gridspec_kw={"height_ratios": [2, 1]})

    # Store the label str expressions
    # micromolar_conc_to_math_exp(vit_conc1, 0) + " original"
    label_line1 = f"{vit_conc1:.1e} $\mu M$ standard"
    label_line2 = f"{vit_conc2:.1e} $\mu M$ standard"
    label_line3 = f"{vit_conc1:.1e} $\mu M$ changing"
    label_line4 = f"{vit_conc2:.1e} $\mu M$ changing"

    # First normal product lines
    line1, = ax1.plot(time1[::plot_di], p1_standard[::plot_di],
                      label=label_line1, color="#F3758A")
    line2, = ax1.plot(time2[::plot_di], p2_standard[::plot_di],
                      label=label_line2, color="#4FD590")
    line3, = ax1.plot(time1[::plot_di], p1[0, ::plot_di],
                      label=label_line3, color="#9B0138")
    line4, = ax1.plot(time2[::plot_di], p2[0, ::plot_di],
                      label=label_line4, color="#057D54")

    fill1 = ax1.fill_between(
        time1[::plot_di], p1_standard[::plot_di], p2_standard[::plot_di], color="#FFCE3A", alpha=0.25)
    fill2 = ax1.fill_between(
        time1[::plot_di], p1[0, ::plot_di], p2[0, ::plot_di], color="#FFCE3A", alpha=0.75)

    # Area line
    area_line1, = ax2.plot(k_D_all, area_array, color="#F3758A")
    area_line2, = ax2.plot(k_D_all[0], area_array[0], color="#9B0138")
    # Plot the position of the k_D that is the standard area
    ax2.scatter(standard_k_D, 1.0, color="#057D54")

    k_D_math_exp = f"$\mathrm{{K}}_{{\mathrm{{D}}}}$: {standard_k_D:.2e} $\mathrm{{\mu M}}$"
    k_D_text = ax1.text(0.65, 0.1, k_D_math_exp, transform=ax1.transAxes,
                        fontsize=10, bbox=dict(facecolor="#FFCE3A", alpha=0.5, boxstyle="round"))

    def init():
        # ax1.set_title("Product concentration over time")
        # ax1.legend(title="Vitamin concentration")
        ax1.legend()
        ax1.set_xlabel(r"Time $[\mathrm{{s}}}$")
        ax1.set_ylabel(r"Product concentration $[\mathrm{{\mu M}}}$")
        ax1.set_xlim(0, t_tot)

        # ax2.set_title("Area fold change between two graphs over time")
        ax2.set_xlabel(r"Dissociation constant $[\mathrm{{\mu M}}}$")
        ax2.set_ylabel("Area fold change")

        # Determine how to set the x-axis for ax2
        x_lim_diff = np.amax(k_D_all)*0.05
        # ax2.set_xlim(np.amin(k_D_all)-x_lim_diff,
        #              np.amax(k_D_all)+x_lim_diff)
        ax2.set_xscale("log")

        # Determine how to set the y-axis for ax2
        y_lim_diff = np.amax(area_array)*0.05
        ax2.set_ylim(np.amin(area_array)-y_lim_diff,
                     np.amax(area_array)+y_lim_diff)
        ax2.yaxis.set_major_locator(MultipleLocator(0.5))

        # Plot vertical and horizontal lines to the the position of the
        # k_D that determines the standard area
        ax2.vlines(standard_k_D, 0, 1, color="#057D54", linestyle="dashed")
        ax2.axhline(1, color="#057D54", linestyle="dashed")
        # ax2.hlines(1, np.amin(k_D_all)-x_lim_diff, standard_k_D,
        #            color="#057D54", linestyle="dashed")

        # Set the character labels
        ax1.text(-0.05, 1.05, "a", transform=ax1.transAxes,
                 size=16, weight="bold")
        ax2.text(-0.05, 1.05, "b", transform=ax2.transAxes,
                 size=16, weight="bold")

        # Set a tight_layout for the figure. This needs to be done
        # after the axis names and title have been set
        # fig.tight_layout()
        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, k_D_text,

    def update(index):

        k_D_math_exp = f"$\mathrm{{K}}_{{\mathrm{{D}}}}$: {k_D_all[index]:.2e} $\mathrm{{\mu M}}$"
        k_D_text.set_text(k_D_math_exp)

        line3.set_data(time1[::plot_di], p1[index, ::plot_di])
        line4.set_data(time2[::plot_di], p2[index, ::plot_di])

        area_line2.set_data(k_D_all[:index], area_array[:index])

        ax1.collections.clear()
        fill1 = ax1.fill_between(
            time1[::plot_di], p1_standard[::plot_di], p2_standard[::plot_di], color="#FFCE3A", alpha=0.4)
        fill2 = ax1.fill_between(
            time1[::plot_di], p1[index, ::plot_di], p2[index, ::plot_di], color="#FFCE3A", alpha=0.8)

        return line1, line2, line3, line4, fill1, fill2, area_line1, area_line2, k_D_text,

    anim = FuncAnimation(fig, update, frames=np.arange(num_steps),
                         init_func=init, blit=False)

    if save_path is not None:
        writermp4 = FFMpegFileWriter(fps=10, bitrate=5000)
        anim.save(f"{save_path}", writer=writermp4)

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
                     label=label_line1, color="#F3758A")
    line2, = ax.plot(time1, frac_umrna_vit[0, :],
                     label=label_line2, color="#4FD590")
    line3, = ax.plot(time1, frac_cmrna[0, :],
                     label=label_line3, color="#9B0138")
    # line4, = ax.plot(time2, absorbance2[0, :],
    #                  label=label_line4, color="#057D54")

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
    #     0.05, 0.09, s_i=1000, low_dna_conc=0.3*10**-4, standard_dna_conc=3*10**-3, high_dna_conc=6*10**-3, num_steps=20, dt=0.01, t_tot=14400)
    # anim_two_vitamin_conc_differing_k_c(
    #     0.05, 0.09, s_i=1000, dna_conc=3*10**-3, low_k_c=(1/60)/10, standard_k_c=1/60, high_k_c=(1/60)*5, num_steps=20, dt=0.01, t_tot=7200)

    # anim_two_vitamin_conc_differing_dna_conc(
    #     0.05, 0.09, s_i=1000, low_dna_conc=0.3*10**-4, standard_dna_conc=3*10**-3, high_dna_conc=6*10**-3, num_steps=361, dt=0.01, t_tot=10800, save_path="modelling/data/plots/T--TUDelft--Animation_Differing_DNA_Conc.mp4")

    anim_two_vitamin_conc_differing_k_c(
        0.05, 0.09, s_i=1000, dna_conc=3*10**-3, low_k_c=(0.017)/10, standard_k_c=0.017, high_k_c=(0.017)*5, num_steps=361, dt=0.01, t_tot=7200, save_path="modelling/data/plots/T--TUDelft--Animation_Differing_k_c.mp4")
