import numpy as np
import matplotlib.pyplot as plt
from models import model_prokaryotic, model_prokaryotic_all
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions
from plot_helpers import custom_aptavita_colors, custom_aptavita_color_cycler, micromolar_conc_to_math_exp


def plot_prokaryotic_different_vitamin_conc(vit_conc_list: list, dna_conc: float = 3*10**-3, s_i: float = 250, save_path: str = None):
    """Function plots the product over time for different vitamin concentrations.

    Parameters
    ----------
    vit_conc_list: list
        all the different vitamin concentrations in micromolar you want to plot
    dna_conc: float
        initial DNA concentrations in the simulation (default 5*10**-3)
    s_i: float
        initial substrate concentration in the simulation (default 150)
    save_path: str
        path to which the figure should be save if the value is not None (default None)
    """
    # Create the custom color cycler
    custom_cycler = custom_aptavita_color_cycler()

    # Create the figure
    fig, ax = plt.subplots(figsize=(12, 6), dpi=125)

    ax.set_prop_cycle(custom_cycler)

    # Store the parameters and constants which are the same
    # for every simulation
    parameters = standard_parameters_prokaryotic()

    # Loop over all the vitamin concentrations and simulate
    # the change in absorbance and plot the proper line in
    # the figure
    for vit_conc in vit_conc_list:
        # Set the initial conditions
        initial_conditions = np.array([dna_conc, s_i, vit_conc])

        # Run the simulation
        time, p = model_prokaryotic(
            parameters, initial_conditions, dt=0.01, t_tot=4800)

        ax.plot(time, p,
                label=f"{int(vit_conc * 1000)} $\mathrm{{nM}}$")

    ax.legend(title="Vitamin concentration")
    # ax.set_title("Absorbance over time for different vitamin concentrations")
    ax.set_xlabel(r"Time $(\mathrm{s})$")
    ax.set_ylabel(r"Product concentration $(\mathrm{\mu M})$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_area_prokaryotic_different_k_c(vit_conc1, vit_conc2, dna_conc: float = 3*10**-3, s_i: float = 250, dt=0.01, t_tot=4800, save_path: str = None):
    """Function plots the absorbance over time for different vitamin concentrations.

    Parameters
    ----------
    vit_conc1:

    vit_conc2:
        all the different vitamin concentrations in micromolar you want to plot
    dna_conc: float
        initial DNA concentrations in the simulation (default 5*10**-3)
    s_i: float
        initial substrate concentration in the simulation (default 150)
    save_path: str
        path to which the figure should be save if the value is not None (default None)
    """
    # Create the custom color cycler
    custom_cycler = custom_aptavita_color_cycler()

    # Create the figure
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(14, 10), dpi=125)

    ax1.set_prop_cycle(custom_cycler)
    ax2.set_prop_cycle(custom_cycler)

    # Store the parameters and constants which are the same
    # for every simulation
    parameters = standard_parameters_prokaryotic()

    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])

    time1, p1 = model_prokaryotic(
        parameters, initial_conditions1, dt=dt, t_tot=t_tot)
    time2, p2 = model_prokaryotic(
        parameters, initial_conditions2, dt=dt, t_tot=t_tot)
    ax1.plot(time1, p1, label=f"{int(vit_conc1 * 1000)} $\mathrm{{nM}}$")
    ax1.plot(time2, p2, label=f"{int(vit_conc2 * 1000)} $\mathrm{{nM}}$")
    ax1.fill_between(time1, p1, p2, color="#FFCE3A", alpha=0.75)

    parameters_adjusted = parameters.copy()
    parameters_adjusted[13] = 1/600

    time_adjusted1, p_adjusted1 = model_prokaryotic(
        parameters_adjusted, initial_conditions1, dt=dt, t_tot=t_tot)
    time_adjusted2, p_adjusted2 = model_prokaryotic(
        parameters_adjusted, initial_conditions2, dt=dt, t_tot=t_tot)
    ax2.plot(time_adjusted1, p_adjusted1,
             label=f"{int(vit_conc1 * 1000)} $\mathrm{{nM}}$")
    ax2.plot(time_adjusted2, p_adjusted2,
             label=f"{int(vit_conc2 * 1000)} $\mathrm{{nM}}$")
    ax2.fill_between(time_adjusted1, p_adjusted1,
                     p_adjusted2, color="#FFCE3A", alpha=0.75)

    # Set ax1 proporties
    ax1.legend(title="Vitamin concentration")
    ax1.set_xlabel(r"Time $(\mathrm{s})$")
    ax1.set_ylabel(r"Product concentration $(\mathrm{\mu M})$")

    # Set ax2 proporties
    ax2.legend(title="Vitamin concentration")
    ax2.set_xlabel(r"Time $(\mathrm{s})$")
    ax2.set_ylabel(r"Product concentration $(\mathrm{\mu M})$")

    # Set the character labels
    ax1.text(-0.05, 1.05, "a", transform=ax1.transAxes,
             size=16, weight="bold")
    ax2.text(-0.05, 1.05, "b", transform=ax2.transAxes,
             size=16, weight="bold")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    # fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_cmrna_diff_vitamin_concentrations_prokaryotic(vit_conc_list: list, dna_conc: float = 5*10**-3, s_i: float = 150, save_path: str = None):
    """Function plots the absorbance over time for different vitamin concentrations.

    Parameters
    ----------
    vit_conc_list: list
        all the different vitamin concentrations in micromolar you want to plot
    dna_conc: float
        initial DNA concentrations in the simulation (default 5*10**-3)
    s_i: float
        initial substrate concentration in the simulation (default 150)
    save_path: str
        path to which the figure should be save if the value is not None (default None)
    """

    # Create the figure
    fig, ax = plt.subplots()

    # Store the parameters and constants which are the same
    # for every simulation
    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()

    # Loop over all the vitamin concentrations and simulate
    # the change in absorbance and plot the proper line in
    # the figure
    for vit_conc in vit_conc_list:
        # Set the initial conditions
        initial_conditions = standard_initial_conditions(
            dna_conc=dna_conc, s_i=s_i, vit_conc=vit_conc)

        # Run the simulation
        results = model_prokaryotic_all(
            parameters, constants, initial_conditions, dt=0.01, t_tot=7200)

        time = results[0]
        cmrna = results[3]

        # Create a label string for the line, this is
        # a string with the correct math expression
        line_label = micromolar_conc_to_math_exp(vit_conc, 2)
        ax.plot(time, cmrna, label=line_label)

    ax.legend(title="Vitamin concentration")
    ax.set_title("Absorbance over time for different vitamin concentrations")
    ax.set_xlabel(r"Time $[s]$")
    ax.set_ylabel(r"Concentration $[\mu M]$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_fraction_mrna_prokaryotic(vit_conc: list, dna_conc: float = 5*10**-3, s_i: float = 150, save_path: str = None):
    # Create the figure
    fig, ax = plt.subplots()

    # Store the parameters and constants which are the same
    # for every simulation
    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()
    initial_conditions = standard_initial_conditions(
        dna_conc=dna_conc, s_i=s_i, vit_conc=vit_conc)
    # Run the simulation
    results = model_prokaryotic_all(
        parameters, constants, initial_conditions, dt=0.01, t_tot=7200)

    time = results[0]

    umrna = results[2]
    umrna_vit = results[6]
    cmrna = results[3]
    total = umrna + umrna_vit + cmrna

    frac_umrna = umrna / total
    frac_umrna_vit = umrna_vit / total
    frac_cmrna = cmrna / total

    ax.plot(time, frac_umrna, label="umRNA")
    ax.plot(time, frac_umrna_vit, label="umRNA_Vit")
    ax.plot(time, frac_cmrna, label="cmRNA")

    ax.legend()

    vit_conc_expr = micromolar_conc_to_math_exp(vit_conc, 2)
    title_text = "Fraction for different mRNA's over time for a vitamin concentration of " + vit_conc_expr

    ax.set_title(title_text)
    ax.set_xlabel(r"Time $[s]$")
    ax.set_ylabel(r"Relative abundance $[AU]$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_absolute_mrna_prokaryotic(vit_conc: list, dna_conc: float = 5*10**-3, s_i: float = 150, save_path: str = None):
    # Create the figure
    fig, ax = plt.subplots()

    # Store the parameters and constants which are the same
    # for every simulation
    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()
    initial_conditions = standard_initial_conditions(
        dna_conc=dna_conc, s_i=s_i, vit_conc=vit_conc)
    # Run the simulation
    results = model_prokaryotic_all(
        parameters, constants, initial_conditions, dt=0.01, t_tot=7200)

    time = results[0]

    umrna = results[2]
    umrna_vit = results[6]
    cmrna = results[3]

    ax.plot(time, umrna, label="umRNA")
    ax.plot(time, umrna_vit, label="umRNA_Vit")
    ax.plot(time, cmrna, label="cmRNA")

    ax.legend()

    vit_conc_expr = micromolar_conc_to_math_exp(vit_conc, 2)
    title_text = "Concentration of different mRNA's over time for a vitamin concentration of " + vit_conc_expr

    ax.set_title(title_text)
    ax.set_xlabel(r"Time $[s]$")
    ax.set_ylabel(r"Concentration $[\mu M]$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_total_absolute_umrna_prokaryotic_differing_k_c(k_c_list: list, k_c_list_names: list, vit_conc: float, dna_conc: float = 3*10**-3, s_i: float = 250, save_path: str = None):

    custom_cycler = custom_aptavita_color_cycler()

    # Create the figure
    fig, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    ax1.set_prop_cycle(custom_cycler)

    for i in range(len(k_c_list)):
        parameters = standard_parameters_prokaryotic()
        parameters[13] = k_c_list[i]
        initial_conditions = np.array([dna_conc, s_i, vit_conc])

        # Run the simulation
        results = model_prokaryotic_all(
            parameters, initial_conditions, dt=0.01, t_tot=7200)

        time = results[0]

        umrna = results[2]
        umrna_vit = results[6]
        total = umrna + umrna_vit

        ax1.plot(time, total*1000,
                 label=f"{k_c_list_names[i]} $\mathrm{{1/s}}$")

    ax1.legend(title="Cleaving rate $(k_{\mathrm{c}})$")

    ax1.set_xlabel(r"Time $(\mathrm{s})$")
    ax1.set_ylabel(r"Concentration umRNA $(\mathrm{{nM}})$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_total_and_bound_absolute_umrna_prokaryotic_differing_k_c(k_c_list: list, k_c_list_names: list, vit_conc: float, dna_conc: float = 3*10**-3, s_i: float = 250, save_path: str = None):

    custom_colors = custom_aptavita_colors()

    # Create the figure
    fig, ax1 = plt.subplots(nrows=1, figsize=(14, 8), dpi=125)
    # ax1.set_prop_cycle(custom_cycler)
    # ax2.set_prop_cycle(custom_cycler)

    for i in range(len(k_c_list)):
        parameters = standard_parameters_prokaryotic()
        parameters[13] = k_c_list[i]
        initial_conditions = np.array([dna_conc, s_i, vit_conc])

        # Run the simulation
        results = model_prokaryotic_all(
            parameters, initial_conditions, dt=0.01, t_tot=7200)

        time = results[0]

        umrna = results[2]
        umrna_vit = results[6]
        total = umrna + umrna_vit

        ax1.plot(time, total*1000,
                 label=f"$\mathrm{{umRNA}}$ $+$ $\mathrm{{umRNA}} \cdot \mathrm{{Vit}}$ $(k_\mathrm{{c}}=$ {k_c_list_names[i]} $\mathrm{{1/s}})$", color=custom_colors[i])
        ax1.plot(time, umrna_vit*1000,
                 label=f"$\mathrm{{umRNA}} \cdot \mathrm{{Vit}}$ $(k_\mathrm{{c}}=$ {k_c_list_names[i]} $\mathrm{{1/s}})$", color=custom_colors[i], alpha=0.5)

    ax1.legend()  # title="Cleaving rate $(k_{\mathrm{c}})$")

    ax1.set_xlabel(r"Time $(\mathrm{s})$")
    ax1.set_ylabel(r"Concentration $(\mathrm{{nM}})$")

    # # Setting proporties for ax2
    # ax2.legend(title="Cleaving rate $(k_{\mathrm{c}})$")
    # ax2.set_xlabel(r"Time $(\mathrm{s})$")
    # ax2.set_ylabel(
    #     r"$\mathrm{{umRNA}} \cdot \mathrm{{Vit}}$ $(\mathrm{{nM}})$")
    # ax2.set_ylim(ax1_ylim)

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


def plot_dumrna_dt_different_dna_conc(dna_conc_list: list, dna_conc_list_names: list, vit_conc_list: float, s_i: float = 250, save_path: str = None):

    custom_cycler = custom_aptavita_color_cycler()
    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    ax1.set_prop_cycle(custom_cycler)

    dt = 0.01
    t_tot = 1800

    parameters = standard_parameters_prokaryotic()

    for j in range(len(vit_conc_list)):
        vit_conc = vit_conc_list[j]
        for i in range(len(dna_conc_list)):
            initial_conditions = np.array([dna_conc_list[i], s_i, vit_conc])
            results = model_prokaryotic_all(
                parameters, initial_conditions, dt=dt, t_tot=t_tot)
            time = results[0]
            cmrna = results[3]
            cmrna_temp = np.concatenate((np.array([0]), cmrna[:-1]))
            dcmrna_dt = (cmrna - cmrna_temp) / dt

            ax1.plot(time, dcmrna_dt*1000,
                     label=fr"cmRNA formation rate $(\mathrm{{DNA}} = {dna_conc_list_names[i]}$ $\mathrm{{nM}})$ $int({vit_conc} *1000)$")

            production_dt = results[-1] / dt
            ax1.plot(time, production_dt*1000, label="umRNA production rate")

    ax1.set_xlabel("Time $(\mathrm{s})$")
    ax1.set_ylabel("$\mathrm{d}C / \mathrm{d}t$ $(\mathrm{nM/s})$")
    ax1.legend()
    plt.show()


def plot_dumrna_dt_different_k_c(k_c_list: list, k_c_list_names: list, vit_conc: float, dna_conc: float = 3*10**-3, s_i: float = 250, save_path: str = None):

    custom_cycler = custom_aptavita_color_cycler()
    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    ax1.set_prop_cycle(custom_cycler)

    dt = 0.01
    t_tot = 7200

    parameters = standard_parameters_prokaryotic()

    initial_conditions = np.array([dna_conc, s_i, vit_conc])

    for i in range(len(k_c_list)):
        parameters[13] = k_c_list[i]
        results = model_prokaryotic_all(
            parameters, initial_conditions, dt=dt, t_tot=t_tot)
        time = results[0]
        cmrna = results[3]
        cmrna_temp = np.concatenate((np.array([0]), cmrna[:-1]))
        dcmrna_dt = (cmrna - cmrna_temp) / dt

        ax1.plot(time, dcmrna_dt*1000,
                 label=f"$\mathrm{{d}}\mathrm{{cmRNA}}/\mathrm{{d}}t$ $(k_{{\mathrm{{c}}}} = {k_c_list_names[i]}$ $1/\mathrm{{s}})$")

        if i == len(k_c_list)-1:
            production_dt = results[-1] / dt
            ax1.plot(time, production_dt*1000,
                     label=f"$\mathrm{{d}}\mathrm{{umRNA}}/\mathrm{{d}}t$", color="#4D94EF")

    ax1.set_xlabel("Time $(\mathrm{s})$")
    ax1.set_ylabel("$\mathrm{d}C / \mathrm{d}t$ $(\mathrm{nM/s})$")
    ax1.legend()
    fig1.tight_layout()
    plt.show()


def k_D_plot():
    vit_tot1 = 0.03
    vit_tot2 = 0.05
    vit_tot3 = 0.07
    vit_tot4 = 0.09
    k_D = 0.05
    umrna = np.linspace(0, 0.02, num=100, dtype=np.float64)
    umrna_vit1 = (vit_tot1*umrna)/(k_D + umrna)
    umrna_vit2 = (vit_tot2*umrna)/(k_D + umrna)
    umrna_vit3 = (vit_tot3*umrna)/(k_D + umrna)
    umrna_vit4 = (vit_tot4*umrna)/(k_D + umrna)

    test = np.linspace(0, 0.022, num=100, dtype=np.float64)

    fig1, ax1 = plt.subplots()
    ax1.plot(umrna, umrna_vit1, label="30")
    ax1.plot(umrna, umrna_vit2, label="50")
    ax1.plot(umrna, umrna_vit3, label="70")
    ax1.plot(umrna, umrna_vit4, label="90")
    ax1.plot(test, test, label="linear line for reference")
    ax1.set_xlim(0, 0.022)
    ax1.set_ylim(0, 0.022)
    ax1.legend()
    plt.show()


if __name__ == "__main__":
    pass
    # plot_prokaryotic_different_vitamin_conc(
    #     [0.05, 0.09], dna_conc=3*10**-3, s_i=250)  # , save_path="test_plot.svg")
    # plot_area_prokaryotic_different_k_c(
    #     0.05, 0.09, dna_conc=3*10**-3, s_i=250, save_path=None)
    # plot_fraction_mrna_prokaryotic(2, save_path="test_plot.svg")
    # plot_total_absolute_mrna_prokaryotic(2, save_path="test_plot2.svg")
    # plot_total_and_bound_absolute_umrna_prokaryotic_differing_k_c(k_c_list=[1/60, 1/300, 1/600], k_c_list_names=[
    #                                                               "1/60", "1/300", "1/600"], vit_conc=0.07, dna_conc=3*10**-3, s_i=250, save_path=None)
    # plot_all_different_absolute_mrna_prokaryotic(2, save_path="test_plot2.svg")
    # plot_cmrna_diff_vitamin_concentrations_prokaryotic(np.linspace(0, 20, 11), save_path="test_plot3.svg")

    # plot_absolute_umrna_prokaryotic_differing_vit_conc(
    #     [1/60, 1/300, 1/600], ["1/60", "1/300", "1/600"], vit_conc=0.05, dna_conc=3*10**-3, s_i=250)

    plot_dumrna_dt_different_k_c(k_c_list=[1/60, 1/300, 1/600], k_c_list_names=[
                                 "1/60", "1/300", "1/600"], vit_conc=0.07, dna_conc=3*10**-3, s_i=250, save_path=None)
    # plot_dumrna_dt_different_dna_conc(dna_conc_list=[
    #                                   3*10**-3, 1*10**-3, 0.3*10**-3], dna_conc_list_names=[
    #     "3", "1", "0.3"], vit_conc_list=[0.05, 0.11], s_i=250, save_path=None)
    # plot_area_prokaryotic_different_k_c(
    #     0.05, 0.09, dna_conc=3*10**-3, s_i=250, dt=0.01, t_tot=4800)
