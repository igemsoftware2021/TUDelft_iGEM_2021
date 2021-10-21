import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from models_area import model_prokaryotic_area
from models import model_prokaryotic, model_prokaryotic_all
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions
from plot_helpers import custom_aptavita_colors, custom_aptavita_color_cycler, micromolar_conc_to_math_exp


def yfp_plot(save_path: str = None):

    # Gain 100
    custom_colors = custom_aptavita_colors()
    time = np.array([0, 300, 600, 900, 1200, 1500, 1800, 2100.1, 2400, 2700.1, 3000.1, 3300.1, 3600.1, 3900.1, 4200.1, 4500.1, 4800.1, 5100.1, 5400.1,
                    5700.1])
    rlu = np.array([117, 593, 2764, 7011, 12349, 18037, 23534, 28223, 32296, 35357, 37888, 39927, 41354, 42476, 43344, 43846, 44502, 45048,
                   45088, 45288])

    time_h = time/3600

    fig1, ax1 = plt.subplots(figsize=(10, 6), dpi=150)
    ax1.plot(time_h, rlu, color=custom_colors[2], linewidth=2)
    # Set minor and major tick locator
    ax1.xaxis.set_major_locator(MultipleLocator(0.25))
    ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax1.set_xlim(0, 1.5)
    ax1.set_xlabel(r"Time $[\mathrm{h}]$")
    ax1.set_ylabel(r"RLU")

    ax1_xlim = ax1.get_xlim()
    start_steady = 0.25  # start time of the steady state regime
    end_steady = 0.6     # end time of the steady state regime

    # Colour the area for every regime
    ax1.axvspan(0, start_steady, color=custom_colors[0], alpha=0.15)
    ax1.axvspan(start_steady, end_steady, color=custom_colors[0], alpha=0.1)
    ax1.axvspan(end_steady, ax1_xlim[1], color=custom_colors[0], alpha=0.05)

    # Add black dashed lines at positions x
    ax1.axvline(start_steady, color="#000000", linestyle="dashed", linewidth=1)
    ax1.axvline(end_steady, color="#000000", linestyle="dashed", linewidth=1)

    # Set the character labels
    ax1.text(0.0555, 0.92, "i", transform=ax1.transAxes,
             size=16, weight="bold")
    ax1.text(0.222, 0.92, "ii", transform=ax1.transAxes,
             size=16, weight="bold")
    ax1.text(0.455, 0.92, "iii", transform=ax1.transAxes,
             size=16, weight="bold")

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def yfp_plot_purefrex2(save_path: str = None):

    custom_colors = custom_aptavita_colors()
    time = np.array([0, 300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000, 3300, 3600, 3900, 4200.1, 4500.1, 4800.1, 5100.1, 5400.1, 5700.1, 6000.1, 6300.1, 6600.1, 6900.1, 7200.1, 7500.1, 7800.1, 8100.1, 8400.1, 8700.1, 9000.1, 9300.1, 9600.1, 9900.1, 10200, 10500, 10800, 11100, 11400, 11700, 12000, 12300, 12600, 12900, 13200, 13500, 13800, 14100, 14400, 14700, 15000, 15300, 15600, 15900, 16200, 16500, 16800, 17100, 17400, 17700, 18000, 18300, 18600, 18900, 19200, 19500, 19800, 20100, 20400, 20700, 21000, 21300, 21600, 21900, 22200, 22500, 22800, 23100, 23400, 23700, 24000, 24300, 24600, 24900, 25200, 25500, 25800, 26100, 26400, 26700, 27000, 27300, 27600, 27900, 28200, 28500,
                    28800, 29101, 29400, 29701, 30001, 30301, 30601, 30901, 31201, 31501, 31801, 32101, 32401, 32701, 33001, 33301, 33601, 33901, 34201, 34501, 34801, 35101, 35401, 35701, 36001, 36301, 36601, 36901, 37201, 37501, 37801, 38101, 38401, 38701, 39001, 39301, 39601, 39901, 40201, 40501, 40801, 41101, 41401, 41701, 42001, 42301, 42601, 42901, 43201, 43501, 43801, 44101, 44401, 44701, 45001, 45301, 45601, 45901, 46201, 46501, 46801, 47101, 47401, 47701, 48001, 48301, 48601, 48901, 49201, 49501, 49801, 50101, 50401, 50701, 51001, 51301, 51601, 51901, 52201, 52501, 52801, 53101, 53401, 53701, 54001, 54301, 54601, 54901, 55201, 55501, 55801, 56101, 56401, 56701, 57001, 57301, 57601])
    rlu = np.array([261, 799, 1673, 2770, 3935, 5076, 6156, 7154, 8108, 8934, 9727, 10412, 11032, 11580, 12087, 12539, 12907, 13296, 13609, 13912, 14185, 14449, 14704, 14944, 15163, 15352, 15524, 15663, 15820, 15994, 16139, 16235, 16373, 16498, 16626, 16747, 16825, 16882, 16988, 17036, 17145, 17180, 17276, 17280, 17338, 17409, 17414, 17480, 17460, 17540, 17541, 17594, 17647, 17682, 17744, 17735, 17776, 17790, 17822, 17798, 17861, 17854, 17921, 17900, 17905, 17954, 17954, 17972, 17973, 18008, 18022, 18000, 18034, 18015, 18058, 18047, 18052, 18060, 18026, 17997, 17986, 18012, 18001, 17975, 17982, 17966, 17965, 17954, 17947, 17939, 17973, 17948, 17960, 17968, 17937, 17958, 17932,
                   17945, 17995, 17970, 18003, 17945, 17965, 17967, 17983, 17957, 17975, 17960, 17983, 17968, 17973, 17979, 17929, 17950, 17959, 17975, 17995, 18000, 18023, 18022, 17988, 18008, 18003, 18023, 17976, 18063, 18036, 18014, 18017, 17979, 18054, 18062, 18016, 18033, 18039, 18039, 18071, 18029, 18050, 18028, 18035, 18040, 18042, 18021, 18039, 18061, 18046, 18040, 18014, 18023, 18026, 18038, 18035, 18020, 18045, 18042, 18014, 18019, 18017, 18010, 18023, 18001, 18007, 17989, 17972, 17970, 17984, 17958, 17987, 17965, 17876, 17929, 17949, 17946, 17900, 17895, 17875, 17930, 17906, 17908, 17922, 17899, 17939, 17997, 17979, 17955, 17970, 17995, 17967, 18018, 17998, 18017, 18011])

    time_h = time/3600

    fig1, ax1 = plt.subplots(figsize=(10, 6), dpi=150)
    ax1.plot(time_h, rlu, color=custom_colors[2], linewidth=2)
    # Set minor and major tick locator
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax1.set_xlim(0, 6)
    ax1.set_xlabel(r"Time $[\mathrm{h}]$")
    ax1.set_ylabel(r"RLU")

    start_steady = 0.17  # start time of the steady state regime
    end_steady = 0.8     # end time of the steady state regime

    # Colour the area for every regime
    ax1.axvspan(0, start_steady, color=custom_colors[0], alpha=0.15)
    ax1.axvspan(start_steady, end_steady, color=custom_colors[0], alpha=0.1)
    ax1.axvspan(end_steady, 6, color=custom_colors[0], alpha=0.05)

    # Add black dashed lines at positions x
    ax1.axvline(start_steady, color="#000000", linestyle="dashed", linewidth=1)
    ax1.axvline(end_steady, color="#000000", linestyle="dashed", linewidth=1)

    # Set the character labels
    ax1.text(0.0075, 0.92, "i", transform=ax1.transAxes,
             size=16, weight="bold")
    ax1.text(0.05, 0.92, "ii", transform=ax1.transAxes,
             size=16, weight="bold")
    ax1.text(0.165, 0.92, "iii", transform=ax1.transAxes,
             size=16, weight="bold")

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_prokaryotic_different_vitamin_conc(vit_conc_list: list, dna_conc: float = 3*10**-3, s_i: float = 1000, dt=0.01, t_tot=4800, save_path: str = None):
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
    fig1, ax1 = plt.subplots(figsize=(10, 6), dpi=150)

    ax1.set_prop_cycle(custom_cycler)

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
            parameters, initial_conditions, dt=dt, t_tot=t_tot)

        ax1.plot(time, p,
                 label=f"{int(vit_conc * 1000)} $\mathrm{{nM}}$")

    ax1.legend(title="Vitamin concentration")
    # ax.set_title("Absorbance over time for different vitamin concentrations")
    ax1.set_xlabel(r"Time $[\mathrm{s}]$")
    ax1.set_ylabel(r"Product concentration $[\mathrm{\mu M}]$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig1.tight_layout()

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_TlR_over_time(vit_conc: float = 0.07, dna_conc: float = 5*10**-3, s_i: float = 1000, save_path: str = None):

    custom_cycler = custom_aptavita_color_cycler()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)

    ax1.set_prop_cycle(custom_cycler)

    parameters = standard_parameters_prokaryotic()

    initial_conditions = np.array([dna_conc, s_i, vit_conc])

    results = model_prokaryotic_all(
        parameters, initial_conditions, dt=0.01, t_tot=43200)

    time = results[0]
    tlr = results[8]

    time_h = time / 3600
    ax1.plot(time_h, tlr, linewidth=2)

    # Set minor and major tick locator
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax1.set_xlim(0, 12)
    ax1.set_xlabel(r"Time $[\mathrm{h}]$")
    # ax1.set_ylabel()
    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_enzyme_mon_and_enzyme_conc(dna_conc: float = 5*10**-3, s_i: float = 1000, save_path: str = None):

    custom_cycler = custom_aptavita_color_cycler()

    fig1, ax1 = plt.subplots(figsize=(10, 6), dpi=150)

    ax1.set_prop_cycle(custom_cycler)

    parameters = standard_parameters_prokaryotic()
    # Set the cleaving rate extremely high so that all the
    # umRNA will be immediately cleaved
    parameters[13] = 1

    # Set vit_conc to 0, since you want complete cleavage of all the umRNA
    vit_conc = 0
    initial_conditions = np.array([dna_conc, s_i, vit_conc])

    results = model_prokaryotic_all(
        parameters, initial_conditions, dt=0.01, t_tot=43200)

    time = results[0]
    enzyme_mon_conc = results[9]
    enzyme_conc = results[10]

    time_h = time / 3600
    ax1.plot(time_h, enzyme_mon_conc, linewidth=2,
             label="$\mathrm{E}^{\\ast}$")
    ax1.plot(time_h, enzyme_conc, linewidth=2, label="$\mathrm{E}$")
    # Set minor and major tick locator
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax1.set_xlim(0, 12)
    ax1.set_xlabel(r"Time $[\mathrm{h}]$")
    ax1.set_ylabel(r"Concentration $[\mathrm{\mu M}]$")
    ax1.legend(title="Species")
    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_area_prokaryotic_different_k_c(vit_conc1, vit_conc2, dna_conc: float = 3*10**-3, s_i: float = 1000, dt=0.01, t_tot=4800, save_path: str = None):
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
    # Plot step, was needed otherwise svg file becomes too large
    plot_step = 100

    # Create the custom color cycler
    custom_cycler = custom_aptavita_color_cycler()

    # Create the figure
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), dpi=150)

    ax1.set_prop_cycle(custom_cycler)
    ax2.set_prop_cycle(custom_cycler)

    # Store the parameters and constants which are the same
    # for every simulation
    parameters = standard_parameters_prokaryotic()

    area_reference = model_prokaryotic_area(
        parameters, dna_conc=dna_conc, s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot)

    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])
    time1, p1 = model_prokaryotic(
        parameters, initial_conditions1, dt=dt, t_tot=t_tot)
    time2, p2 = model_prokaryotic(
        parameters, initial_conditions2, dt=dt, t_tot=t_tot)
    ax1.plot(time1[::plot_step], p1[::plot_step],
             label=f"{int(vit_conc1 * 1000)} $\mathrm{{nM}}$")
    ax1.plot(time2[::plot_step], p2[::plot_step],
             label=f"{int(vit_conc2 * 1000)} $\mathrm{{nM}}$")
    ax1.fill_between(time1[::plot_step], p1[::plot_step],
                     p2[::plot_step], color="#FFCE3A", alpha=0.75)

    parameters_adjusted = parameters.copy()
    parameters_adjusted[13] = 1/600

    area_adjusted = model_prokaryotic_area(
        parameters_adjusted, dna_conc=dna_conc, s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot)

    fold_change = (area_reference + (area_adjusted -
                   area_reference)) / area_reference
    print(f"Fold change: {fold_change}")

    time_adjusted1, p_adjusted1 = model_prokaryotic(
        parameters_adjusted, initial_conditions1, dt=dt, t_tot=t_tot)
    time_adjusted2, p_adjusted2 = model_prokaryotic(
        parameters_adjusted, initial_conditions2, dt=dt, t_tot=t_tot)
    ax2.plot(time_adjusted1[::plot_step], p_adjusted1[::plot_step],
             label=f"{int(vit_conc1 * 1000)} $\mathrm{{nM}}$")
    ax2.plot(time_adjusted2[::plot_step], p_adjusted2[::plot_step],
             label=f"{int(vit_conc2 * 1000)} $\mathrm{{nM}}$")
    ax2.fill_between(time_adjusted1[::plot_step], p_adjusted1[::plot_step],
                     p_adjusted2[::plot_step], color="#FFCE3A", alpha=0.75)

    # Set ax1 proporties
    ax1.legend(title="Vitamin concentration")
    ax1.set_xlabel(r"Time $[\mathrm{s}]$")
    ax1.set_ylabel(r"Product concentration $[\mathrm{\mu M}]$")

    # Set ax2 proporties
    ax2.legend(title="Vitamin concentration")
    ax2.set_xlabel(r"Time $[\mathrm{s}]$")
    ax2.set_ylabel(r"Product concentration $[\mathrm{\mu M}]$")

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
    else:
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
    else:
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
    else:
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
    else:
        plt.show()


def plot_total_absolute_umrna_prokaryotic_differing_k_c(k_c_list: list, k_c_list_names: list, vit_conc: float, dna_conc: float = 3*10**-3, s_i: float = 1000, save_path: str = None):

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
                 label=f"{k_c_list_names[i]} $\mathrm{{s}}^{{-1}}$")

    ax1.legend(title="Cleaving rate $(k_{\mathrm{c}})$")

    ax1.set_xlabel(r"Time $[\mathrm{s}]$")
    ax1.set_ylabel(r"Concentration umRNA $[\mathrm{{nM}}]$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_total_and_bound_absolute_umrna_prokaryotic_differing_k_c(k_c_list: list, k_c_list_names: list, vit_conc: float, dna_conc: float = 3*10**-3, s_i: float = 1000, save_path: str = None):

    custom_colors = custom_aptavita_colors()

    # Create the figure
    fig, ax1 = plt.subplots(nrows=1, figsize=(14, 8), dpi=150)
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

        ax1.plot(time[::10], total[::10]*1000,
                 label=f"$\mathrm{{umRNA}} + \mathrm{{umRNA}} \cdot \mathrm{{Vit}}\;(k_\mathrm{{c}}= {k_c_list_names[i]}\;\mathrm{{s}}^{{-1}})$", color=custom_colors[i])
        ax1.plot(time[::10], umrna[::10]*1000,
                 label=f"$\mathrm{{umRNA}}\;(k_\mathrm{{c}} = {k_c_list_names[i]}\;\mathrm{{s}}^{{-1}})$", color=custom_colors[i], alpha=0.5)

    ax1.legend()

    ax1.set_xlabel(r"Time $[\mathrm{s}]$")
    ax1.set_ylabel(r"Concentration $[\mathrm{{nM}}]$")

    # # Setting proporties for ax2
    # ax2.legend(title="Cleaving rate $(k_{\mathrm{c}})$")
    # ax2.set_xlabel(r"Time $[\mathrm{s}]$")
    # ax2.set_ylabel(
    #     r"$\mathrm{{umRNA}} \cdot \mathrm{{Vit}}\;[\mathrm{{nM}}]$")
    # ax2.set_ylim(ax1_ylim)

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_dumrna_dt_different_dna_conc(dna_conc_list: list, dna_conc_list_names: list, vit_conc_list: float, s_i: float = 1000, save_path: str = None):

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
                     label=fr"cmRNA formation rate $(\mathrm{{DNA}} = {dna_conc_list_names[i]}\;\mathrm{{nM}})\;int({vit_conc} *1000)$")

            production_dt = results[-1] / dt
            ax1.plot(time, production_dt*1000, label="umRNA production rate")

    ax1.set_xlabel("Time $[\mathrm{s}]$")
    ax1.set_ylabel("$\mathrm{d}C / \mathrm{d}t\;[\mathrm{nM/s}]$")
    ax1.legend()
    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_dumrna_dt_different_k_c(k_c_list: list, k_c_list_names: list, vit_conc: float, dna_conc: float = 3*10**-3, s_i: float = 1000, save_path: str = None):

    custom_colors = custom_aptavita_colors()

    fig1, ax1 = plt.subplots(figsize=(12, 7), dpi=150)

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

        ax1.plot(time, dcmrna_dt*1000, color=custom_colors[i], alpha=0.8, zorder=1,
                 label=f"$\mathrm{{d}}[\mathrm{{cmRNA}}]/\mathrm{{d}}t\;(k_{{\mathrm{{c}}}} = {k_c_list_names[i]}\;\mathrm{{s}}^{{-1}})$")

        dproduction_dt = results[-1] / dt

        if i == len(k_c_list)-1:
            dproduction_dt = results[-1] / dt
            ax1.plot(time, dproduction_dt*1000, zorder=1,
                     label=f"$\mathrm{{d}}[\mathrm{{umRNA}}]/\mathrm{{d}}t$", color="#4D94EF", alpha=0.8)

        difference = dproduction_dt - dcmrna_dt
        difference[:50] = 1000
        min_idx = np.argmin(difference)
        ax1.scatter(time[min_idx], dcmrna_dt[min_idx] * 1000, marker="o",
                    color=custom_colors[i], s=30, alpha=1, zorder=2, edgecolors="#000000")

    ax1.set_xlabel("Time $[\mathrm{s}]$")
    ax1.set_ylabel("$\mathrm{d}C / \mathrm{d}t\;[\mathrm{nM/s}]$")
    ax1.legend()
    # ax1.tight_layout()

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
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

    # plot_prokaryotic_different_vitamin_conc(
    #     [0.07], dna_conc=3*10**-3, s_i=1000, dt=0.01, t_tot=5400, save_path="modelling/data/plots/T--TUDelft--Model_Example_Plot.svg")

    # plot_prokaryotic_different_vitamin_conc(
    #     [0.05, 0.09], dna_conc=3*10**-3, s_i=1000, dt=0.01, t_tot=5400, save_path="modelling/data/plots/T--TUDelft--Model_Two_Vitamin_Plot.svg")

    # plot_area_prokaryotic_different_k_c(
    #     0.05, 0.09, dna_conc=3*10**-3, s_i=1000, dt=0.01, t_tot=7200, save_path="modelling/data/plots/T--TUDelft--Area_Example_Plot.svg")

    # plot_dumrna_dt_different_k_c(k_c_list=[0.017, 0.0033, 0.0017], k_c_list_names=[
    #                              "0.017", "0.0033", "0.0017"], vit_conc=0.07, dna_conc=3*10**-3, s_i=1000, save_path="modelling/data/plots/T--TUDelft--dumRNA_dt_Different_k_c_Plot.svg")

    plot_total_and_bound_absolute_umrna_prokaryotic_differing_k_c(k_c_list=[0.017, 0.0033, 0.0017], k_c_list_names=[
                                                                  "0.017", "0.0033", "0.0017"], vit_conc=0.07, dna_conc=3*10**-3, s_i=1000, save_path="modelling/data/plots/T--TUDelft--umRNA_Bound_Unbound_Different_k_c_Plot.svg")

    # yfp_plot(save_path="modelling/data/plots/T--TUDelft--YFP_Plot_Three_Regimes.svg")
    # plot_enzyme_mon_and_enzyme_conc(
    #     dna_conc=3*10**-3, s_i=1000)  # , save_path="modelling/data/plots/T--TUDelft--Plot_Enzyme_Concentration.svg")
    # plot_TlR_over_time(vit_conc=0.07, dna_conc=5 *
    #                    10**-3, s_i=1000, save_path=None)
