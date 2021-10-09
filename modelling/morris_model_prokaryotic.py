import time
import numpy as np
from morris_method import morris_analysis_prokaryotic, morris_problem_description_prokaryotic
from morris_method import morris_analysis_area_prokaryotic, morris_area_problem_description_prokaryotic
from morris_method import morris_datawriter
from standard_values import standard_parameters_prokaryotic


def run_morris_analysis_prokaryotic(path="modelling/data", tag=f"_{int(time.time())}"):
    # Set Morris method values
    num_trajectories = 1000
    num_levels = 4
    num_optimal_trajectories = 70

    # Set the time of the simulation and timesteps
    t_tot = 10800
    dt = 0.01

    var_names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
                 "k_tlr", "k_m", "deg_mrna", "deg_tlr", "k_on", "k_off", "k_c", "dna_conc", "vit_conc"]

    # Determine the number of parameters
    num_parameters = len(var_names)

    # Set the parameters and the parameter ranges
    parameters = standard_parameters_prokaryotic()

    dna_conc = 0.3*10**-3
    s_i = 250
    vit_conc = 0.07

    initial_conditions = np.array([dna_conc, s_i, vit_conc])

    lower_range = 1/3
    upper_range = 3

    # Problem definition for prokaryotic system
    prokaryotic_problem = {
        "num_vars": num_parameters,
        "names": var_names,
        "bounds": [[lower_range * parameters[0], upper_range * parameters[0]],     # (0) k_ts
                   [lower_range * parameters[1], upper_range * \
                    parameters[1]],     # (1) k_tl
                   [lower_range * parameters[2], upper_range * \
                    parameters[2]],     # (2) k_mat
                   [lower_range * parameters[3], upper_range * \
                    parameters[3]],     # (3) k_cat
                   [lower_range * parameters[4], upper_range * \
                    parameters[4]],     # (4) k_s
                   [lower_range * parameters[5], upper_range * \
                    parameters[5]],     # (5) kc_s
                   [lower_range * parameters[6], upper_range * \
                    parameters[6]],     # (6) k_l
                   [lower_range * parameters[7], upper_range * \
                    parameters[7]],     # (7) k_tlr
                   [lower_range * parameters[8], upper_range * \
                    parameters[8]],     # (8) k_m
                   [lower_range * parameters[9], upper_range * \
                    parameters[9]],     # (9) deg_mrna
                   [lower_range * parameters[10], upper_range * \
                    parameters[10]],  # (10) deg_tlr
                   [lower_range * parameters[11], upper_range * \
                    parameters[11]],  # (11) k_on
                   [lower_range * parameters[12], upper_range * \
                    parameters[12]],  # (12) k_off
                   [lower_range * parameters[13], upper_range * \
                    parameters[13]],  # ,  # (13) k_c
                   [lower_range * dna_conc, upper_range * \
                    dna_conc],  # (14) dna_conc
                   [lower_range * vit_conc, upper_range * vit_conc]]  # (15) vit_conc
    }

    # Check if the number of variable names givens is equal to the number of bounds given
    if len(prokaryotic_problem["names"]) != len(prokaryotic_problem["bounds"]):
        raise ValueError(
            "The number of variables is not equal to the number of bounds")

    # Write a description
    morris_problem_description_prokaryotic(
        s_i, prokaryotic_problem, num_trajectories, num_levels, num_optimal_trajectories, dt=dt, t_tot=t_tot, path=path, tag=tag)

    (time, mu, mu_star, sigma, mu_star_conf_level) = morris_analysis_prokaryotic(
        prokaryotic_problem, num_trajectories, initial_conditions, dt=dt, t_tot=t_tot, num_levels=num_levels, optimal_trajectories=num_optimal_trajectories)

    morris_datawriter(prokaryotic_problem, mu, mu_star, sigma,
                      mu_star_conf_level, time=time, path=path, tag=tag)


def run_morris_analysis_area_prokaryotic(path="modelling/data", tag=f"_{int(time.time())}"):
    # Set Morris method values
    num_trajectories = 1000
    num_levels = 4
    num_optimal_trajectories = 70

    # Set the time of the simulation and timesteps
    t_tot = 7200
    dt = 0.01

    var_names = ["k_ts", "k_tl", "k_mat", "k_cat", "k_s", "kc_s", "k_l",
                 "k_tlr", "k_m", "deg_mrna", "deg_tlr", "k_on", "k_off", "k_c", "dna_conc"]  # , "vit_conc"]

    # Determine the number of parameters
    num_parameters = len(var_names)

    # Set the parameters and the parameter ranges
    parameters = standard_parameters_prokaryotic()

    dna_conc = 3*10**-3
    s_i = 250
    vit_conc1 = 0.05
    vit_conc2 = 0.09

    lower_range = 1/3
    upper_range = 3

    # Problem definition for prokaryotic system
    prokaryotic_problem = {
        "num_vars": num_parameters,
        "names": var_names,
        "bounds": [[lower_range * parameters[0], upper_range * parameters[0]],     # (0) k_ts
                   [lower_range * parameters[1], upper_range * \
                    parameters[1]],     # (1) k_tl
                   [lower_range * parameters[2], upper_range * \
                    parameters[2]],     # (2) k_mat
                   [lower_range * parameters[3], upper_range * \
                    parameters[3]],     # (3) k_cat
                   [lower_range * parameters[4], upper_range * \
                    parameters[4]],     # (4) k_s
                   [lower_range * parameters[5], upper_range * \
                    parameters[5]],     # (5) kc_s
                   [lower_range * parameters[6], upper_range * \
                    parameters[6]],     # (6) k_l
                   [lower_range * parameters[7], upper_range * \
                    parameters[7]],     # (7) k_tlr
                   [lower_range * parameters[8], upper_range * \
                    parameters[8]],     # (8) k_m
                   [lower_range * parameters[9], upper_range * \
                    parameters[9]],     # (9) deg_mrna
                   [lower_range * parameters[10], upper_range * \
                    parameters[10]],  # (10) deg_tlr
                   [lower_range * parameters[11], upper_range * \
                    parameters[11]],  # (11) k_on
                   [lower_range * parameters[12], upper_range * \
                    parameters[12]],  # (12) k_off
                   [lower_range * parameters[13], upper_range * \
                    parameters[13]],  # ,  # (13) k_c
                   [lower_range * dna_conc, upper_range * \
                    dna_conc]],  # (14) dna_conc
    }

    # Check if the number of variable names givens is equal to the number of bounds given
    if len(prokaryotic_problem["names"]) != len(prokaryotic_problem["bounds"]):
        raise ValueError(
            "The number of variables is not equal to the number of bounds")

    # Write a description
    morris_area_problem_description_prokaryotic(
        s_i, vit_conc1, vit_conc2, prokaryotic_problem, num_trajectories, num_levels, num_optimal_trajectories, dt=dt, t_tot=t_tot, path=path, tag=tag)

    (mu, mu_star, sigma, mu_star_conf_level) = morris_analysis_area_prokaryotic(prokaryotic_problem,
                                                                                num_trajectories, dna_conc=dna_conc, s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot, num_levels=num_levels, optimal_trajectories=num_optimal_trajectories)

    morris_datawriter(prokaryotic_problem, mu, mu_star, sigma,
                      mu_star_conf_level, path=path, tag=tag)


if __name__ == "__main__":
    pass
    # run_morris_analysis_prokaryotic()
    run_morris_analysis_area_prokaryotic()
