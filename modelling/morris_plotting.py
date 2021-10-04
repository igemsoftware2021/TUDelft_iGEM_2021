import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from morris_method import morris_datareader
from plot_helpers import custom_aptavita_colors, custom_aptavita_color_cycler, add_value_labels


def plot_morris_analysis(path="modelling/data", tag="_1633190385", save_path=None):
    parameters, data_dict = morris_datareader(path=path, tag=tag)

    custom_cycler = custom_aptavita_color_cycler()

    fig1, ax1 = plt.subplots(figsize=(20, 12), dpi=100)
    fig2, ax2 = plt.subplots(figsize=(20, 12), dpi=100)
    fig3, ax3 = plt.subplots(figsize=(20, 12), dpi=100)
    fig4, ax4 = plt.subplots(figsize=(20, 12), dpi=100)

    # Set the color cycler
    ax1.set_prop_cycle(custom_cycler)
    ax2.set_prop_cycle(custom_cycler)
    ax3.set_prop_cycle(custom_cycler)
    ax4.set_prop_cycle(custom_cycler)

    for i in range(len(parameters)):
        ax1.plot(data_dict["time"][::10], data_dict["mu"]
                 [::10, i], label=parameters[i])
        # ax2.plot(data_dict["time"][::10], data_dict["mu"]
        #          [::10, i], label=parameters[i])
        # ax3.plot(data_dict["time"][::10], data_dict["mu"]
        #          [::10, i], label=parameters[i])
        ax2.plot(data_dict["time"][::10], data_dict["mu_star"]
                 [::10, i], label=parameters[i])
        ax3.plot(data_dict["time"][::10], data_dict["sigma"]
                 [::10, i], label=parameters[i])
        ax4.plot(data_dict["time"][::10], data_dict["mu_star_conf_level"]
                 [::10, i], label=parameters[i], alpha=0.5)

    # Set all proporties for ax1
    ax1.legend()
    ax1.set_xlabel(r"Time $\mathrm{(s)}$")
    ax1.set_ylabel(r"$\mathrm{\mu}$ $(\mathrm{\mu M})$")

    # Set all proporties for ax2
    ax2.legend()
    ax2.set_xlabel(r"Time $\mathrm{(s)}$")
    ax2.set_ylabel(r"$\mathrm{\mu}^{\ast}$ $(\mathrm{\mu M})$")

    # Set all proporites for ax3
    ax3.legend()
    ax3.set_xlabel(r"Time $\mathrm{(s)}$")
    ax3.set_ylabel(r"$\mathrm{\sigma}$ $(\mathrm{\mu M})$")

    ax4.legend()
    ax4.set_xlabel(r"Time $\mathrm{(s)}$")
    ax4.set_ylabel(
        r"$\mathrm{\mu^{\ast}}$ $95\%$-confidence interval $(\mathrm{\mu M})$")

    if save_path is not None:
        fig1.savefig(f"{save_path}/plot_mu{tag}", format="svg", dpi=1200)
        fig2.savefig(f"{save_path}/plot_mu_star{tag}", format="svg", dpi=1200)
        fig3.savefig(f"{save_path}/plot_sigma{tag}", format="svg", dpi=1200)
        fig4.savefig(f"{save_path}/plot_mu_star_ci{tag}",
                     format="svg", dpi=1200)

    plt.show()


def plot_morris_analysis_area(path="modelling/data", tag="_1633293118", save_path=None):
    parameters, data_dict = morris_datareader(path=path, tag=tag)

    mu = data_dict["mu"].reshape(-1)
    mu_star = data_dict["mu_star"].reshape(-1)
    sigma = data_dict["sigma"].reshape(-1)

    custom_colors = custom_aptavita_colors()

    fig1, ax1 = plt.subplots(figsize=(16, 9), dpi=125)
    fig2, ax2 = plt.subplots(figsize=(16, 9), dpi=125)
    fig3, ax3 = plt.subplots(figsize=(16, 9), dpi=125)

    ax1.bar(np.arange(0, mu.shape[0]), mu, color=custom_colors)
    ax2.bar(np.arange(0, mu_star.shape[0]), mu_star, color=custom_colors)
    ax3.bar(np.arange(0, sigma.shape[0]), sigma, color=custom_colors)

    # Set all proporties for ax1
    # ax1.set_xlabel(r"Time $(s)$")
    ax1.set_xticks(np.arange(0, len(parameters)))
    ax1.set_xticklabels(parameters)
    ax1.set_ylabel(r"$\mathrm{\mu}$")
    ax1.yaxis.set_major_locator(MultipleLocator(5000))

    # Set all proporties for ax2
    ax2.set_xticks(np.arange(0, len(parameters)))
    ax2.set_xticklabels(parameters)
    ax2.set_ylabel(r"$\mathrm{\mu}^{\ast}$")
    ax2.yaxis.set_major_locator(MultipleLocator(5000))

    # Set all proporites for ax3
    ax3.set_xticks(np.arange(0, len(parameters)))
    ax3.set_xticklabels(parameters)
    ax3.set_ylabel(r"$\mathrm{\sigma}$")
    ax3.yaxis.set_major_locator(MultipleLocator(5000))

    if save_path is not None:
        fig1.savefig(f"{save_path}/plot_mu{tag}", format="svg", dpi=1200)
        fig2.savefig(f"{save_path}/plot_mu_star{tag}", format="svg", dpi=1200)
        fig3.savefig(f"{save_path}/plot_sigma{tag}", format="svg", dpi=1200)

    plt.show()


if __name__ == "__main__":
    # plot_morris_analysis(path="modelling/data",
    #                      tag="_1633190385")  # , save_path="modelling/data/plots")
    plot_morris_analysis_area(
        path="modelling/data")  # , tag="_1633293196")  # , save_path="modelling/data/plots")
