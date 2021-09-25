import matplotlib.pyplot as plt
from models import model_prokaryotic_readout
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions
from plot_helpers import micromolar_conc_to_math_exp


def plot_vitamin_concentrations_prokaryotic(vit_conc_list: list, dna_conc: float = 5*10**-3, s_i: float = 150, save_path: str = None):
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
        time, absorbance = model_prokaryotic_readout(
            parameters, constants, initial_conditions, dt=0.01, t_tot=7200)

        # Create a label string for the line, this is
        # a string with the correct math expression
        line_label = micromolar_conc_to_math_exp(vit_conc, 2)
        ax.plot(time, absorbance, label=line_label)

    ax.legend()
    ax.set_title("Absorbance over time for different vitamin concentrations")
    ax.set_xlabel(r"Time $[s]$")
    ax.set_ylabel(r"Absorbance $[AU]$")

    # Tight layout should be set after plotting and setting the
    # correct title and labels
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(f"{save_path}", format="svg", dpi=1200)

    plt.show()


if __name__ == "__main__":
    plot_vitamin_concentrations_prokaryotic(
        [5, 10, 20], save_path="test_plot.svg")
