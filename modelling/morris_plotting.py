import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from morris_method import morris_datareader
from plot_helpers import custom_aptavita_colors, custom_aptavita_color_cycler, senstivity_analysis_factor_names
from standard_values import standard_parameters_prokaryotic
from models_area import model_prokaryotic_area


def plot_morris_analysis(path="modelling/data", tag="_1633190385", save_path=None):
    parameters, data_dict = morris_datareader(path=path, tag=tag)

    custom_cycler = custom_aptavita_color_cycler()
    factor_names = senstivity_analysis_factor_names()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    fig2, ax2 = plt.subplots(figsize=(14, 8), dpi=125)
    fig3, ax3 = plt.subplots(figsize=(14, 8), dpi=125)
    fig4, ax4 = plt.subplots(figsize=(14, 8), dpi=125)

    # Set the color cycler
    ax1.set_prop_cycle(custom_cycler)
    ax2.set_prop_cycle(custom_cycler)
    ax3.set_prop_cycle(custom_cycler)
    ax4.set_prop_cycle(custom_cycler)

    for i in range(len(parameters)):
        ax1.plot(data_dict["time"][::10], data_dict["mu"]
                 [::10, i], label=factor_names[i])
        # ax2.plot(data_dict["time"][::10], data_dict["mu"]
        #          [::10, i], label=parameters[i])
        # ax3.plot(data_dict["time"][::10], data_dict["mu"]
        #          [::10, i], label=parameters[i])
        ax2.plot(data_dict["time"][::10], data_dict["mu_star"]
                 [::10, i], label=factor_names[i])
        ax3.plot(data_dict["time"][::10], data_dict["sigma"]
                 [::10, i], label=factor_names[i])
        ax4.plot(data_dict["time"][::10], data_dict["mu_star_conf_level"]
                 [::10, i], label=factor_names[i], alpha=0.5)

    # Set all proporties for ax1
    ax1.legend()
    ax1.set_xlabel(r"Time $\mathrm{(s)}$")
    ax1.set_ylabel(r"$\mathrm{{\mu}}$ $[\mathrm{\mu M}]$")

    # Set all proporties for ax2
    ax2.legend()
    ax2.set_xlabel(r"Time $\mathrm{(s)}$")
    ax2.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M}]$")

    # Set all proporites for ax3
    ax3.legend()
    ax3.set_xlabel(r"Time $\mathrm{(s)}$")
    ax3.set_ylabel(r"$\mathrm{{\sigma}}$ $[\mathrm{\mu M}]$")

    ax4.legend()
    ax4.set_xlabel(r"Time $\mathrm{(s)}$")
    ax4.set_ylabel(
        r"$\mathrm{{\mu}}^{\ast}$ $95\%$-confidence interval $[\mathrm{\mu M}]$")

    if save_path is not None:
        fig1.savefig(f"{save_path}/plot_mu{tag}", format="svg", dpi=1200)
        fig2.savefig(f"{save_path}/plot_mu_star{tag}", format="svg", dpi=1200)
        fig3.savefig(f"{save_path}/plot_sigma{tag}", format="svg", dpi=1200)
        fig4.savefig(f"{save_path}/plot_mu_star_ci{tag}",
                     format="svg", dpi=1200)

    plt.show()


def plot_morris_analysis_mu_star_subplots(path="modelling/data", tag="_1633190385", save_path=None):

    fill_plot_step = 10
    plot_steps = 10
    plot_time = 72001

    parameters, data_dict = morris_datareader(
        path=path, tag=tag, data_names=["mu_star", "mu_star_conf"])

    custom_cycler = custom_aptavita_color_cycler()
    custom_colors = custom_aptavita_colors()
    factor_names = senstivity_analysis_factor_names()

    fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows=2, ncols=2, figsize=(14, 10), dpi=125)

    # Set the color cycler
    # ax1.set_prop_cycle(custom_cycler)
    ax2.set_prop_cycle(custom_cycler)
    ax3.set_prop_cycle(custom_cycler)
    ax4.set_prop_cycle(custom_cycler)

    for i in range(len(parameters)):
        x = data_dict["time"][:plot_time:plot_steps]
        y = data_dict["mu_star"][:plot_time:plot_steps, i]
        ci = data_dict["mu_star_conf"][:plot_time:plot_steps, i]

        if parameters[i] in {"deg_mrna", "kc_s", "k_tlr"}:
            ax1.plot(
                x, y, color=custom_colors[0], label="Insensitive" if parameters[i] == "deg_mrna" else "")
            ax1.fill_between(x[::fill_plot_step], (y-ci)[::fill_plot_step], (y+ci)[::fill_plot_step],
                             color=custom_colors[0], alpha=0.05)

            ax2.plot(x, y, label=factor_names[i])

            # Plot the category as color also in the ax3, ax4
            ax3.plot(x, y, color="#755F26", alpha=0.25)
            ax4.plot(x, y, color="#755F26", alpha=0.25)

        if parameters[i] in {"deg_tlr", "k_m", "k_on", "k_off", "k_c", "vit_conc"}:

            ax1.plot(
                x, y, color=custom_colors[1], label="Sensitive" if parameters[i] == "k_m" else "")
            ax1.fill_between(x[::fill_plot_step], (y-ci)[::fill_plot_step], (y+ci)[::fill_plot_step],
                             color=custom_colors[1], alpha=0.05)

            ax3.plot(x, y, label=factor_names[i])

            # Plot the category as color also in the ax2, ax4
            ax2.plot(x, y, color="#755F26", alpha=0.25)
            ax4.plot(x, y, color="#755F26", alpha=0.25)

        if parameters[i] in {"dna_conc", "k_ts", "k_s", "k_tl", "k_l", "k_mat", "k_cat"}:
            ax1.plot(
                x, y, color=custom_colors[2], label="Very sensitive" if parameters[i] == "k_ts" else "")
            ax1.fill_between(x[::fill_plot_step], (y-ci)[::fill_plot_step], (y+ci)[::fill_plot_step],
                             color=custom_colors[2], alpha=0.05)

            ax4.plot(x, y, label=factor_names[i])

            # Plot the category as color also in the ax2, ax3
            ax2.plot(x, y, color="#755F26", alpha=0.25)
            ax3.plot(x, y, color="#755F26", alpha=0.25)

    # Needed to add the label "other" at the end of the legend
    ax2.plot([], [], color="#755F26", alpha=0.25, label="Other")
    ax3.plot([], [], color="#755F26", alpha=0.25, label="Other")
    ax4.plot([], [], color="#755F26", alpha=0.25, label="Other")

    # Set all proporties for ax1
    ax1.legend(loc="upper left")
    ax1.set_xlabel(r"Time $[\mathrm{s}]$")
    ax1.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M}]$")
    ax1_ylim = ax1.get_ylim()

    # Set all proporties for ax2
    ax2.legend()
    ax2.set_xlabel(r"Time $[\mathrm{s}]$")
    ax2.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M}]$")
    # ax2.set_ylim(ax1_ylim)

    # Set all proporites for ax3
    ax3.legend()
    ax3.set_xlabel(r"Time $[\mathrm{s}]$")
    ax3.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M}]$")
    # ax3.set_ylim(ax1_ylim)

    ax4.legend()
    ax4.set_xlabel(r"Time $[\mathrm{s}]$")
    ax4.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M}]$")
    # ax4.set_ylim(ax1_ylim)

    # Set the character labels
    ax1.text(-0.05, 1.05, "a", transform=ax1.transAxes,
             size=16, weight="bold")
    ax2.text(-0.05, 1.05, "b", transform=ax2.transAxes,
             size=16, weight="bold")
    ax3.text(-0.05, 1.05, "c", transform=ax3.transAxes,
             size=16, weight="bold")
    ax4.text(-0.05, 1.05, "d", transform=ax4.transAxes,
             size=16, weight="bold")

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_morris_analysis_area(path="modelling/data", tag="_1633293118", save_path=None):
    parameters, data_dict = morris_datareader(path=path, tag=tag)

    mu = data_dict["mu"].reshape(-1)
    mu_star = data_dict["mu_star"].reshape(-1)
    sigma = data_dict["sigma"].reshape(-1)

    custom_colors = custom_aptavita_colors()
    factor_names = senstivity_analysis_factor_names()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    fig2, ax2 = plt.subplots(figsize=(14, 8), dpi=125)
    fig3, ax3 = plt.subplots(figsize=(14, 8), dpi=125)

    ax1.bar(np.arange(0, mu.shape[0]), mu, color=custom_colors)
    ax2.bar(np.arange(0, mu_star.shape[0]), mu_star, color=custom_colors)
    ax3.bar(np.arange(0, sigma.shape[0]), sigma, color=custom_colors)

    # Set all proporties for ax1
    # ax1.set_xlabel(r"Time $[s]$")
    ax1.set_xticks(np.arange(0, len(parameters)))
    ax1.set_xticklabels(factor_names[:len(parameters)])
    ax1.set_ylabel(r"$\mathrm{{\mu}}$ $[\mathrm{\mu M \cdot s}]$")
    ax1.yaxis.set_major_locator(MultipleLocator(5000))

    # Set all proporties for ax2
    ax2.set_xticks(np.arange(0, len(parameters)))
    ax2.set_xticklabels(factor_names[:len(parameters)])
    ax2.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M \cdot s}]$")
    ax2.yaxis.set_major_locator(MultipleLocator(5000))

    # Set all proporites for ax3
    ax3.set_xticks(np.arange(0, len(parameters)))
    ax3.set_xticklabels(factor_names[:len(parameters)])
    ax3.set_ylabel(r"$\mathrm{{\sigma}}$ $[\mathrm{\mu M \cdot s}]$")
    ax3.yaxis.set_major_locator(MultipleLocator(5000))

    if save_path is not None:
        fig1.savefig(f"{save_path}_mu.svg", format="svg", dpi=1200)
        fig2.savefig(f"{save_path}_mu_star.svg", format="svg", dpi=1200)
        fig3.savefig(f"{save_path}_sigma.svg", format="svg", dpi=1200)
    else:
        plt.show()


def plot_morris_analysis_area_fold_change(path="modelling/data", tag="_1633293118", save_path=None, standard_area=8758.052481272032):
    parameters, data_dict = morris_datareader(path=path, tag=tag)

    mu = data_dict["mu"].reshape(-1)
    mu_star = data_dict["mu_star"].reshape(-1)

    custom_colors = custom_aptavita_colors()
    factor_names = senstivity_analysis_factor_names()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    fig2, ax2 = plt.subplots(figsize=(14, 8), dpi=125)
    fig3, (ax3, ax4) = plt.subplots(
        nrows=2, ncols=1, figsize=(12, 10), dpi=150)

    mu_fc = (standard_area + mu)/standard_area
    mu_star_fc = (standard_area + mu_star)/standard_area

    # The mu fold change will be offset by 1, since that is the standard
    # Hence all the values need to be lowered by 1 to get the correct
    # bars
    mu_fc_bar_value = mu_fc - 1

    ax1.bar(np.arange(0, mu_fc.shape[0]),
            mu_fc_bar_value, bottom=1, color=custom_colors)
    ax2.bar(np.arange(0, mu_star_fc.shape[0]), mu_star_fc, color=custom_colors)

    ax3.bar(np.arange(0, mu_fc.shape[0]),
            mu_fc_bar_value, bottom=1, color=custom_colors)
    ax4.bar(np.arange(0, mu_star_fc.shape[0]), mu_star_fc, color=custom_colors)

    # Set all proporties for ax1
    ax1.set_xticks(np.arange(0, len(parameters)))
    ax1.set_xticklabels(factor_names[:len(parameters)])
    ax1.set_ylabel(r"Fold change")
    ax1.yaxis.set_minor_locator(MultipleLocator(0.1))

    # Set all proporties for ax2
    ax2_ylim = ax2.get_ylim()
    ax2.set_ylim(1, ax2_ylim[1])
    ax2.set_xticks(np.arange(0, len(parameters)))
    ax2.set_xticklabels(factor_names[:len(parameters)])
    ax2.set_ylabel(
        r"Fold change")
    ax2.yaxis.set_minor_locator(MultipleLocator(0.1))

    # Set all proporties for ax3
    ax3.set_xticks(np.arange(0, len(parameters)))
    ax3.set_xticklabels(factor_names[:len(parameters)])
    ax3.set_ylabel(r"Fold change")
    ax3.yaxis.set_minor_locator(MultipleLocator(0.1))

    # Set all proporties for ax4
    ax4_ylim = ax2.get_ylim()
    ax4.set_ylim(1, ax2_ylim[1])
    ax4.set_xticks(np.arange(0, len(parameters)))
    ax4.set_xticklabels(factor_names[:len(parameters)])
    ax4.set_ylabel(
        r"Fold change")
    ax4.yaxis.set_minor_locator(MultipleLocator(0.1))

    # Set the character labels
    ax3.text(-0.05, 1.05, "a", transform=ax3.transAxes,
             size=16, weight="bold")
    ax4.text(-0.05, 1.05, "b", transform=ax4.transAxes,
             size=16, weight="bold")

    if save_path is not None:
        fig1.savefig(f"{save_path}_Mu_FC{tag}.svg", format="svg", dpi=1200)
        fig2.savefig(f"{save_path}_Mu_Star_FC{tag}.svg",
                     format="svg", dpi=1200)
        fig3.savefig(f"{save_path}_Mu_FC_Mu_Star_FC_Subplot{tag}.svg",
                     format="svg", dpi=1200)
    else:
        plt.show()


def plot_morris_analysis_area_fold_change_multiple_DNA_conc(save_path=None):
    # Ran for 3 hours
    # parameters1, data_dict1 = morris_datareader(
    #     path="modelling/data", tag="_1634025739")
    # parameters2, data_dict2 = morris_datareader(
    #     path="modelling/data", tag="_1634030488")
    # parameters3, data_dict3 = morris_datareader(
    #     path="modelling/data", tag="_1634030607")
    # parameters4, data_dict4 = morris_datareader(
    #     path="modelling/data", tag="_1634030685")

    dt = 0.01
    t_tot = 18000

    # Ran for 5 hours
    parameters1, data_dict1 = morris_datareader(
        path="modelling/data", tag="_1634031996")  # 3 nM
    parameters2, data_dict2 = morris_datareader(
        path="modelling/data", tag="_1634032760")  # 1 nM
    parameters3, data_dict3 = morris_datareader(
        path="modelling/data", tag="_1634032582")  # 0.3 nM
    parameters4, data_dict4 = morris_datareader(
        path="modelling/data", tag="_1634033312")  # 0.1 nM

    mu1 = data_dict1["mu"].reshape(-1)
    mu2 = data_dict2["mu"].reshape(-1)
    mu3 = data_dict3["mu"].reshape(-1)
    mu4 = data_dict4["mu"].reshape(-1)

    standard_parameter_values = standard_parameters_prokaryotic()
    standard_area1 = model_prokaryotic_area(
        standard_parameter_values, 3*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)
    standard_area2 = model_prokaryotic_area(
        standard_parameter_values, 1*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)
    standard_area3 = model_prokaryotic_area(
        standard_parameter_values, 0.3*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)
    standard_area4 = model_prokaryotic_area(
        standard_parameter_values, 0.1*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)

    mu_fc1 = mu1 / standard_area1
    mu_fc2 = mu2 / standard_area2
    mu_fc3 = mu3 / standard_area3
    mu_fc4 = mu4 / standard_area4

    # Not used
    # mu_star = data_dict["mu_star"].reshape(-1)
    # sigma = data_dict["sigma"].reshape(-1)

    custom_colors = custom_aptavita_colors()
    factor_names = senstivity_analysis_factor_names()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)

    # Position of bars on x-axis
    ind = np.arange(mu_fc1.shape[0])

    # Width of a bar
    width = 0.20

    ax1.bar(ind-1.5*width, mu_fc1, width,
            color=custom_colors[0], label="DNA 3 nM")
    ax1.bar(ind-0.5*width, mu_fc2, width,
            color=custom_colors[1], label="DNA 1 nM")
    ax1.bar(ind+0.5*width, mu_fc3, width,
            color=custom_colors[2], label="DNA 0.3 nM")
    ax1.bar(ind+1.5*width, mu_fc4, width,
            color=custom_colors[3], label="DNA 0.1 nM")
    # ax2.bar(np.arange(0, mu_star_fc.shape[0]), mu_star_fc, color=custom_colors)
    # ax3.bar(np.arange(0, sigma_fc.shape[0]), sigma_fc, color=custom_colors)

    # Set all proporties for ax1
    # ax1.set_xlabel(r"Time $[s]$")
    ax1.legend()
    ax1.set_xticks(np.arange(0, len(parameters1)))
    ax1.set_xticklabels(factor_names[:len(parameters1)])
    ax1.set_ylabel(r"Fold change")
    # ax1.set_ylabel(r"$\mathrm{{\mu}}$ $[\mathrm{\mu M \cdot s}]$")
    ax1.yaxis.set_major_locator(MultipleLocator(0.5))

    if save_path is not None:
        fig1.savefig(f"{save_path}_mu_fc.svg", format="svg", dpi=1200)
    else:
        plt.show()


def plot_morris_analysis_area_mu_star_fold_change_multiple_DNA_conc(save_path=None):
    # Ran for 3 hours
    # parameters1, data_dict1 = morris_datareader(
    #     path="modelling/data", tag="_1634025739")
    # parameters2, data_dict2 = morris_datareader(
    #     path="modelling/data", tag="_1634030488")
    # parameters3, data_dict3 = morris_datareader(
    #     path="modelling/data", tag="_1634030607")
    # parameters4, data_dict4 = morris_datareader(
    #     path="modelling/data", tag="_1634030685")

    dt = 0.01
    t_tot = 18000

    # Ran for 5 hours
    parameters1, data_dict1 = morris_datareader(
        path="modelling/data", tag="_1634031996")  # 3 nM
    parameters2, data_dict2 = morris_datareader(
        path="modelling/data", tag="_1634032760")  # 1 nM
    parameters3, data_dict3 = morris_datareader(
        path="modelling/data", tag="_1634032582")  # 0.3 nM
    parameters4, data_dict4 = morris_datareader(
        path="modelling/data", tag="_1634033312")  # 0.1 nM

    mu_star1 = data_dict1["mu_star"].reshape(-1)
    mu_star2 = data_dict2["mu_star"].reshape(-1)
    mu_star3 = data_dict3["mu_star"].reshape(-1)
    mu_star4 = data_dict4["mu_star"].reshape(-1)

    standard_parameter_values = standard_parameters_prokaryotic()
    standard_area1 = model_prokaryotic_area(
        standard_parameter_values, 3*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)
    standard_area2 = model_prokaryotic_area(
        standard_parameter_values, 1*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)
    standard_area3 = model_prokaryotic_area(
        standard_parameter_values, 0.3*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)
    standard_area4 = model_prokaryotic_area(
        standard_parameter_values, 0.1*10**-3, 250, 0.05, 0.09, dt=dt, t_tot=t_tot)

    mu_star_fc1 = (standard_area1 + mu_star1) / standard_area1
    mu_star_fc2 = (standard_area2 + mu_star2) / standard_area2
    mu_star_fc3 = (standard_area3 + mu_star3) / standard_area3
    mu_star_fc4 = (standard_area4 + mu_star4) / standard_area4

    # Not used
    # mu_star_star = data_dict["mu_star_star"].reshape(-1)
    # sigma = data_dict["sigma"].reshape(-1)

    custom_colors = custom_aptavita_colors()
    factor_names = senstivity_analysis_factor_names()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)

    # Position of bars on x-axis
    ind = np.arange(mu_star_fc1.shape[0])

    # Width of a bar
    width = 0.20

    ax1.bar(ind-1.5*width, mu_star_fc1, width,
            color=custom_colors[0], label="DNA 3 nM")
    ax1.bar(ind-0.5*width, mu_star_fc2, width,
            color=custom_colors[1], label="DNA 1 nM")
    ax1.bar(ind+0.5*width, mu_star_fc3, width,
            color=custom_colors[2], label="DNA 0.3 nM")
    ax1.bar(ind+1.5*width, mu_star_fc4, width,
            color=custom_colors[3], label="DNA 0.1 nM")
    # ax2.bar(np.arange(0, mu_star_star_fc.shape[0]), mu_star_fc, color=custom_colors)
    # ax3.bar(np.arange(0, sigma_fc.shape[0]), sigma_fc, color=custom_colors)

    # Set all proporties for ax1
    # ax1.set_xlabel(r"Time $[s]$")
    ax1.legend()
    ax1.set_xticks(np.arange(0, len(parameters1)))
    ax1.set_xticklabels(factor_names[:len(parameters1)])
    ax1.set_ylabel(r"Fold change")
    # ax1.set_ylabel(r"$\mathrm{{\mu}}$ $[\mathrm{\mu M \cdot s}]$")
    ax1.yaxis.set_major_locator(MultipleLocator(0.5))
    ax1.axhline(1, color="#4D94EF", linestyle="dashed")

    if save_path is not None:
        fig1.savefig(f"{save_path}_mu_star_fc.svg", format="svg", dpi=1200)
    else:
        plt.show()


def morris_method_visualization():
    num_factors = 3
    num_levels = 4
    num_trajectories = 3
    trajectory_length = num_factors + 1

    factor_a = np.arange(num_levels)
    factor_b = np.arange(num_levels)
    factor_c = np.arange(num_levels)

    # Create trajectories
    # Every row is a trajectory step
    # Every column is a seperate factor
    # Every 3d dimension is a full trajectory
    trajectories = np.zeros((trajectory_length, num_factors, num_trajectories))

    # A numpy state that gives a very nice visualization
    np.random.set_state(('MT19937', np.array([2147483648, 3934104647, 1969102952, 4261734483, 2676231918,
                                              1961734023, 3410130820, 2928927076, 1306731865,  789512148,
                                              150769002,  662051233, 2343168163, 2881228086,  343204690,
                                              902510012, 1761393830, 1721288546, 1628181754, 4033512893,
                                              1786262226, 2310580882,  556907703, 2692594442,  837779989,
                                              1974378804, 1839633392, 1090956113, 3409894772,  371300203,
                                              2259790167, 2138388610, 1529118910, 3109964772, 2080982719,
                                              1038539185,   57401473,  706146300, 3093009949, 1348448680,
                                              1726805462,  939096739,  728662136, 1367097771,  205482421,
                                              1494636572, 3621638011, 2350525037,  806264532, 3869362395,
                                              57261087, 3481856095, 2171969119, 4149451513, 2639539794,
                                              1581151049,  938321303,  349267403, 2581687863, 2763364578,
                                              3083783940, 1457017618, 4010148448, 3499546977, 3262895201,
                                              109237176, 1675332589, 3318551450, 3825275499, 2641045249,
                                              1391522102, 3286155305, 1668812482,  535814197, 1630333370,
                                              3642887157, 3238886663, 1085990192, 1844861904, 2537348939,
                                              3590605662, 3517391197, 2634764797,  763200682, 2724115444,
                                              3068802240, 4073567819, 2317930847,  280120392, 3541897351,
                                              3481008685,  414526178, 1016473564, 3855656346, 3415685489,
                                              987349294, 2611813175, 3747546853, 2378620016, 1626511525,
                                              3559360717,  176287339,  960784386, 3432269155, 4121161113,
                                              980338228, 1499481974, 3754929944,  584260697,  326194043,
                                              3646964422, 2604261710,  412413774, 2073145067,  410683916,
                                              1963273938, 4056167918,  580496947, 2736845094, 2816592534,
                                              2548428959, 2884914158, 2561319548, 1257250471, 1070378502,
                                              57624485, 3493209168, 1899588814,  830245478,  837160412,
                                              3124706569, 3220395428, 4250669622, 4172784093, 2992436378,
                                              4169540465,  738413241, 1290560054,  748411661,  307740800,
                                              2887789866, 2334846880, 4111711285,  623991110,   97818510,
                                              2624555925, 1193760977, 2414282288, 2819858105,  190130815,
                                              1587673383, 3010354201,  664475233, 1013419773, 3520407939,
                                              1287383788,  836470161, 2279618730,  277680491, 3168812823,
                                              1766027481, 2007539474, 1437016337,  455441266, 2810472046,
                                              4075401475,  953122953, 4169753982, 4281269780,  870553785,
                                              2292948413, 2642737491,  963300064,  183809534,   89026536,
                                              1683643949, 1966545389,  233089438,  143614457,  736038920,
                                              1483777370, 2352342943, 3308552242, 3212443150, 3255227004,
                                              174015813, 2893503905, 1351235805,  550319623, 2357770050,
                                              3191863641, 3743525590, 1755851954, 2169588331,  625950585,
                                              3034003072, 1134794697, 2497119421, 4111024335,  255140034,
                                              1338713201, 1814508029,  766170940, 1430795167, 1716698200,
                                              2008927287, 2881141660, 1592705689, 2065571908, 4147155402,
                                              2295395878, 1326658244, 2072646688, 1731706653, 3081724201,
                                              1443893607,  407552688, 1993922370,  630915789, 3051959230,
                                              1351711195, 2098375512,  337760455,   83964662, 1651105011,
                                              2242146846,  388535454,  518875834, 3595398311, 2486578908,
                                              2303418893, 3792355325,  698826653,  366552963,  831236509,
                                              2320922241,     121995, 2467415323, 1326686901, 1347264058,
                                              228309057, 1068911820, 2367865172, 3570793218, 1348923980,
                                              3257237647, 1655639721,  167928693,  247425989, 1843985394,
                                              3147603786, 3320737209, 1821164358, 4203003982,  661643529,
                                              632032347, 1889372551, 2709846974, 2794330040, 3074113900,
                                              400579210,  585698691, 2770861537,  626379426, 1080038618,
                                              1844382005, 1968899277, 2908758967,  754394109, 1316173049,
                                              4064469348, 3683585016, 3634661630, 2577411162, 3279491727,
                                              492886199, 1937176540, 1893468150, 1979075696,  129662157,
                                              740439000,  977818322,  850909409,   56047872, 3812372854,
                                              244669407, 3753874521,  505940447, 4021233717, 2452133010,
                                              3422252365,  577705528,  851586344, 3736127322, 3452766091,
                                              3792851941, 3411489901, 1959237021, 2010050776, 3289470150,
                                              294092813,   45799344, 2317926030, 1757355720, 4039535132,
                                              750218957, 3612724259, 2989356538, 3861691654,  346584322,
                                              4123310786, 2750832257, 1289339555, 3438135799,  778620350,
                                              3519311251, 3782690237, 2612585074,  280889816, 3053699060,
                                              2259616949, 3967616327, 3373571362,   94589940, 3904326365,
                                              3501973201, 3037825751, 1807312974, 1017992471, 3118737982,
                                              1038440819, 1949652706,  204731449,  603544179, 3641621288,
                                              2481666972, 2676000474, 3322200830, 1446346435,  924569292,
                                              3498745120, 1136825832, 3360565448, 1233982202, 3229128834,
                                              1093342897,  837488829, 1166759900, 2927225139,  259529522,
                                              229970771, 3174332099, 2024637495, 2214121026,  373599218,
                                              208880676, 4064045782, 1000538953,  480764385, 3728693192,
                                              4141277932, 2185040487,  882982583, 1456670212, 2784743281,
                                              625094984, 1056831835, 1894660761,  474378140, 4269569798,
                                              1982498908,   91545251, 3663473270, 3062756903, 1314173010,
                                              2111455405, 2911239935,  556029246, 3135482933, 3989326239,
                                              945386899, 2249129668, 1849052752, 3195601806, 3647729265,
                                              1118941835, 2410114193, 2974268176, 1292920442, 1645399928,
                                              4276423715, 2430155566,  582316047,  674032374, 1161799337,
                                              2955580157, 1362230331,  340339907, 3428362622, 1725176483,
                                              2035975936, 2348943352, 1628302075, 3157091834, 3917744942,
                                              1642140875, 3996171702, 1024176764, 1555415375, 4192845683,
                                              4037747965,  401813924, 2413129631, 2629570550, 2697968937,
                                              1713980894, 1827563892, 2046895803, 2119872236,  371390069,
                                              2201448157, 1774124628, 3203390121, 2477567579, 1505830139,
                                              3223781244, 3360278879, 2086770870, 2301830841, 3874321895,
                                              506501270, 1884256824, 2946135872, 3688596385,  100323468,
                                              1634826918, 1806914201, 4186613280, 2669641241, 2558691829,
                                              1239303588, 1124245528, 3209392797,  593140493,   16925640,
                                              1137069068, 2954875769,   90980507, 2808768368,  850480523,
                                              2824810089, 2776069920, 4213046135,  552518788, 2740122260,
                                              2759809928, 3568087830, 1409525209, 2342477065, 1947386565,
                                              886584018, 3075561691, 1302028733, 1437428186, 2115053899,
                                              2467777005, 1018291783, 2466340566, 2243173918, 3648582994,
                                              2015955631, 3591487846, 3860144367, 2244915284, 3424903029,
                                              3552736661, 2717927242, 3733091528, 2505467312, 1716933025,
                                              3312648168, 1954881343, 3194966489, 1466237746, 1553181305,
                                              525250980, 2541883588, 2741370514, 4200495204, 1911535541,
                                              3818887049, 3900810326, 3888213835, 2443099028,  998842927,
                                              2035476261, 3900754266, 2120833299,  747799780, 3784068699,
                                              3237121752, 2430535285,  788997889, 4031748408, 1741928105,
                                              3388997658,  775569132, 1038534065, 2751594443, 2397995044,
                                              1878695324,  917950273,  184217568, 4286116162, 1944809181,
                                              4035229395, 1123662739,  752922506, 1372149007, 2633507599,
                                              266952531,  858508273, 1333527274,  597568841,  658828366,
                                              576769258, 1515087208, 1338679229, 1211561068, 3453805606,
                                              1519150970, 4221902450,  177719869, 3086777558, 3747363546,
                                              3476259908,   83107598, 2518813368, 3372552019,  202290777,
                                              4133718057, 2201594483, 1362896602,  142180238, 1354797914,
                                              2210986354, 1886089949, 2683055035, 3308785362, 3038673140,
                                              134577731, 2959516049,   79331200, 2599612328, 3336812021,
                                              2036570548, 1635462385, 3156677933, 1680409743,  836913045,
                                              3301397021, 4269082639, 1721557482, 3584246388, 2336492654,
                                              873365924, 3131064742, 2078009004, 3713443628, 3337748459,
                                              3907540798, 3034953906, 2224317534, 1358306620, 3805730736,
                                              1439534138,  427419491, 2886807851, 2906463499,   43230404,
                                              407029288, 4228292204, 3660460520, 1198232732, 1227562129,
                                              1777684015, 1899045287,  966032195, 3748789115, 2863990591,
                                              540499473, 3195586063, 3822675587, 3533679205,  625127860,
                                              1523262236, 3692045754,  678101569, 2562790471, 3319721641,
                                              1099040120, 2042229283,  361097427,  547570420, 1076380514,
                                              3494852613,  741021389,  289004238, 2462186075, 2532964163,
                                              1247564778,  634277273, 2609602401, 2181279465, 2756328716,
                                              3253377053, 1786112591, 3603493522, 1579989050,  458079130,
                                              3492542113, 3716822400, 3481545762, 2672637242], dtype=np.uint32), 623, 0, 0.0))

    # print(np.random.get_state())

    for i in range(num_trajectories):
        for j in range(trajectory_length):
            if j == 0:
                point = np.array([np.random.choice(factor_a), np.random.choice(
                    factor_b), np.random.choice(factor_c)])
                trajectories[j, :, i] = point
            else:
                prev_point = trajectories[j-1, :, i]
                point = prev_point.copy()
                while True:
                    point = prev_point.copy()
                    random_val = np.random.uniform()
                    if random_val < 1/2:
                        point[j-1] = point[j-1] - 2
                    else:
                        point[j-1] = point[j-1] + 2
                    if (point[j-1] != prev_point[j-1]) and (point[j-1] >= np.amin(factor_a) and point[j-1] <= np.amax(factor_a)):
                        break
                trajectories[j, :, i] = point

    custom_cycler = custom_aptavita_color_cycler()

    fig1 = plt.figure(figsize=(6, 5), dpi=150, frameon=False)
    ax1 = fig1.add_subplot(projection="3d")
    ax1.set_prop_cycle(custom_cycler)
    for i in range(num_trajectories):
        # print(trajectories[:, :, i])
        ax1.plot(trajectories[:, 0, i],
                 trajectories[:, 1, i], trajectories[:, 2, i], linewidth=2)
        ax1.scatter(trajectories[:, 0, i],
                    trajectories[:, 1, i], trajectories[:, 2, i], s=50, alpha=1.0)

    # Plot all x-axis lines
    for i in range(factor_a.shape[0]):
        for k in range(factor_c.shape[0]):
            ax1.plot(np.ones(factor_a.shape) *
                     factor_a[i], factor_b, factor_c[k], color="#000000", linewidth=0.25, alpha=0.5)
            ax1.scatter(np.ones(factor_a.shape) *
                        factor_a[i], factor_b, factor_c[k], color="#000000", s=1, alpha=0.5)

    # Plot all y-axis lines
    for j in range(factor_b.shape[0]):
        for k in range(factor_c.shape[0]):
            ax1.plot(factor_a, np.ones(factor_b.shape) *
                     factor_b[j], factor_c[k], color="#000000", linewidth=0.25, alpha=0.5)
            ax1.scatter(factor_a, np.ones(factor_b.shape) *
                        factor_b[j], factor_c[k], color="#000000", s=1, alpha=0.5)

    # Plot all z-axis lines
    for i in range(factor_a.shape[0]):
        for j in range(factor_b.shape[0]):
            for k in range(factor_c.shape[0]):
                ax1.plot(np.ones(factor_a.shape) * factor_a[i], np.ones(
                    factor_b.shape) * factor_b[j], factor_c, color="#000000", linewidth=0.25, alpha=0.5)
                ax1.scatter(np.ones(factor_a.shape) * factor_a[i], np.ones(
                    factor_b.shape) * factor_b[j], factor_c, color="#000000", s=1, alpha=0.5)

    ax1.grid(False)
    ax1.set_xlim(0, 3)
    ax1.set_ylim(0, 3)
    ax1.set_zlim(0, 3)

    ax1.set_xlabel("Factor A")
    ax1.set_ylabel("Factor B")
    ax1.set_zlabel("Factor C")

    # Turn off ticks
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_zticks([])

    # Transparent spines
    # 0.0, means 0 transparency so it won't be visible
    ax1.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    # 0.0, means 0 transparency so it won't be visible
    ax1.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    # 0.0, means 0 transparency so it won't be visible
    ax1.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Transparent panes
    # 0.0, means 0 transparency so it won't be visible
    ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # 0.0, means 0 transparency so it won't be visible
    ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # 0.0, means 0 transparency so it won't be visible
    ax1.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # This plot should be manually saved, since you need to orient
    # the cube correctly so the trajectories are well visible
    plt.show()


if __name__ == "__main__":
    morris_method_visualization()
    # plot_morris_analysis_mu_star_subplots(
    #     path="modelling/data", tag="_1633520689", save_path="modelling/data/plots/T--TUDelft--Morris_Mu_Star_Subplots_1633520689.svg")
    # plot_morris_analysis_area(
    #     path="modelling/data", tag="_1633591400", save_path="modelling/data/plots/T--TUDelft--Morris_Area_1633591400")

    # parameters = standard_parameters_prokaryotic()
    # standard_area = model_prokaryotic_area(
    #     parameters, 3*10**-3, 250, 0.05, 0.09)
    # plot_morris_analysis_area_fold_change(
    #     path="modelling/data", tag="_1633591400", standard_area=standard_area, save_path="modelling/data/plots/T--TUDelft--Morris_Area")
    # plot_morris_analysis_area_fold_change(
    #     path="modelling/data", tag="_1633591400", save_path=None, standard_area=standard_area)
    # plot_morris_analysis_area_fold_change_multiple_DNA_conc(save_path=None)
    # plot_morris_analysis_area_mu_star_fold_change_multiple_DNA_conc(
    #     save_path=None)
