import csv
import time
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from tqdm import tqdm
import numpy as np


def morris_analysis(problem, trajectories, func, initial_conditions, dt: int = 0.01, t_tot: int = 7200, num_levels: int = 4, optimal_trajectories: int = None, local_optimization: bool = True, num_resamples: int = 1000, conf_level: float = 0.95, print_to_console: bool = False, seed: int = None):
    """Function does Morris analysis on a function/model.

    IMPORTANT
    ---------
    dna_conc and vit_conc are overriden if they are part of the part of the problem!!!


    This is a wrapper function that creates model inputs required for Method of Morris. It then runs these model inputs through
    a model using a parallelized function. Then it analyses the model output and returns the Numpy arrays: time, mu, mu_star, sigma, mu_star_conf_level.

    This is a wrapper function that uses two functions from the SAlib library. Hence the description for the parameters:
    problem, trajectories, num_levels, optimal_trajectories, local_optimization, num_resamples, conf_level, print_to_console
    and seed where directly taken from https://salib.readthedocs.io/en/latest/api.html.

    Parameters
    ----------
    problem: dict
        The problem definition
    trajectories: int
        Number of trajectories to generate.
    func: function
        The function that has as input args: parameters, constants, initial_conditions[, dt, t_tot]
    parameters: numpy.array
        The Numpy array containing all the parameters for the model of dtype=float
    constants: numpy.array
        The Numpy array containing all the constants for the model of dtype=float
    initial_condition: numpy.array
        The Numpy array containing all the initial conditions for the model of dtype=float
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)
    num_levels: int
        The number of grid levels (should be even) (default 4)
    optimal_trajectories: int
        The number of optimal trajectories to sample (between 2 and N) (default None)
    local_optimization: bool
        Flag whether to use local optimization according to Ruano et al. (2012) Speeds up the process
        tremendously for bigger N and num_levels. If set to False brute force method used, unless
        gurobipy is available (default True)
    num_resamples: int
        The number of resamples used to compute the confidence intervals (default 1000)
    conf_level: float
        The confidence interval level (default 0.95)
    print_to_console: bool
        Print results directly to console (default False)
    seed: int
        Seed to generate a random number (default None)

    Returns
    -------
    time: numpy.array
        The Numpy array containing all the timepoints at which the simulation was run. The Numpy array is of dtype=float.
    mu: numpy.array
        The Numpy array containing the mean elementary effects over time. Each row is a timepoint and every column
        contains the mean elementary effects for a certain parameter. The Numpy array is of dtype=float.
    mu_star: numpy.array
        The Numpy array containing the absolute mean elementary effects over time. Each row is a timepoint and every column
        contains the absolute mean elementary effects for a certain parameter. The Numpy array is of dtype=float.
    sigma: numpy.array
        The Numpy array containing the standard deviation sof the mean elementary effect over time. Each row is a timepoint
        and every column contains the standard deviations of the mean elementary effect for a certain parameter.
        The Numpy array is of dtype=float.
    mu_star_conf_level: numpy.array
        The Numpy array containing the bootstrapped confidence intervals. Each row is a timepoint and every column
        contains the bootstrapped confidence interval for a certain parameter. The Numpy array is of dtype=float.
    """
    # First do some checks whether the problem variable contains all the needed keys for the function.
    if "num_vars" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'num_vars' (https://salib.readthedocs.io/en/latest/basics.html)")
    if "names" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'names' (https://salib.readthedocs.io/en/latest/basics.html)")
    if "bounds" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'bounds' (https://salib.readthedocs.io/en/latest/basics.html)")

    # Unpack initial conditions array
    dna_conc, s_i, vit_conc = initial_conditions

    # Determening the timepoints of the simulation.
    n = int(np.ceil(t_tot/dt) + 1)  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints.

    # Retrieve number of parameters
    num_parameters = problem["num_vars"]

    # Defining arrays for sensitivity indices.
    # Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
    mu = np.zeros((n, num_parameters),
                  dtype=np.float32)  # The mean elementary effect
    mu_star = np.zeros((n, num_parameters), dtype=np.float32)
    sigma = np.zeros((n, num_parameters), dtype=np.float32)
    mu_star_conf_level = np.zeros((n, num_parameters), dtype=np.float32)

    # Generating input parameters for the model
    # Each column is a parameter, one row contains all input parameters for one simulation
    model_input = morris_sample.sample(
        problem, trajectories, num_levels=num_levels, optimal_trajectories=optimal_trajectories, local_optimization=local_optimization, seed=seed)

    # Check if you need to override the dna_conc variable then do np.delete()
    # https://stackoverflow.com/questions/24027040/how-to-extract-all-columns-but-one-from-an-array-or-matrix-in-python
    if "dna_conc" in problem["names"]:
        print("dna_conc is part of the problem, so will be overriden!!")
        index = problem["names"].index("dna_conc")
        dna_conc_array = model_input[:, index]
        # Remove the column with the dna concentration, since this is a seperate value from the input parameters
        model_input_temp = np.delete(model_input, index, axis=1)
    else:
        dna_conc_array = np.ones(
            model_input.shape[0], dtype=np.float64) * dna_conc
        model_input_temp = model_input

    # Check if you need to override the vit_conc variable then do np.delete()
    # https://stackoverflow.com/questions/24027040/how-to-extract-all-columns-but-one-from-an-array-or-matrix-in-python
    if "vit_conc" in problem["names"]:
        print("vit_conc is part of the problem, so will be overriden!!")
        index = problem["names"].index("vit_conc")
        vit_conc_array = model_input[:, index]
        # Remove the column with the dna concentration, since this is a seperate value from the input parameters
        model_input_temp = np.delete(model_input, index, axis=1)
    else:
        vit_conc_array = np.ones(
            model_input.shape[0], dtype=np.float64) * vit_conc
        model_input_temp = model_input_temp

    s_i_array = np.ones(
        model_input.shape[0], dtype=np.float64) * s_i

    initial_conditions = np.stack(
        (dna_conc_array, s_i_array, vit_conc_array), axis=1)

    model_output = func(parameters=model_input,
                        initial_conditions=initial_conditions, dt=dt, t_tot=t_tot)

    # Running the Morris analysis at each timepoint (using the output of all the different simulations)
    for ii in tqdm(range(n)):
        indices_dict = morris_analyze.analyze(problem, model_input,
                                              model_output[ii, :], num_levels=num_levels)
        # Each column contains the sensitivity index of one parameter, each column contains the sensitivity indeces at one timestep
        mu[ii, :] = indices_dict['mu']
        mu_star[ii, :] = indices_dict['mu_star']
        sigma[ii, :] = indices_dict['sigma']
        mu_star_conf_level[ii, :] = np.array(
            indices_dict['mu_star_conf'], dtype=np.float32)

    return time, mu, mu_star, sigma, mu_star_conf_level


def morris_analysis_area(problem, trajectories, func, dna_conc, s_i, vit_conc1, vit_conc2, dt: int = 0.01, t_tot: int = 7200, num_levels: int = 4, optimal_trajectories: int = None, local_optimization: bool = True, num_resamples: int = 1000, conf_level: float = 0.95, print_to_console: bool = False, seed: int = None):
    # TODO finish this function, I only copied it and changed the input arguments
    """Function does Morris analysis on a function/model.

    IMPORTANT
    ---------
    dna_conc is overriden if it is part of the problem!!!

    This is a wrapper function that creates model inputs required for Method of Morris. It then runs these model inputs through
    a model using a parallelized function. Then it analyses the model output and returns the Numpy arrays: time, mu, mu_star, sigma, mu_star_conf_level.

    This is a wrapper function that uses two functions from the SAlib library. Hence the description for the parameters:
    problem, trajectories, num_levels, optimal_trajectories, local_optimization, num_resamples, conf_level, print_to_console
    and seed where directly taken from https://salib.readthedocs.io/en/latest/api.html.

    Parameters
    ----------
    problem: dict
        The problem definition
    trajectories: int
        Number of trajectories to generate.
    func: function
        The function that has as input args: parameters, constants, initial_conditions[, dt, t_tot]
    parameters: numpy.array
        The Numpy array containing all the parameters for the model of dtype=float
    constants: numpy.array
        The Numpy array containing all the constants for the model of dtype=float
    initial_condition: numpy.array
        The Numpy array containing all the initial conditions for the model of dtype=float
    dt: int
        The time each timestep takes in seconds. (default 0.01)
    t_tot: int
        The total time the model should run in seconds. (default 7200)
    num_levels: int
        The number of grid levels (should be even) (default 4)
    optimal_trajectories: int
        The number of optimal trajectories to sample (between 2 and N) (default None)
    local_optimization: bool
        Flag whether to use local optimization according to Ruano et al. (2012) Speeds up the process
        tremendously for bigger N and num_levels. If set to False brute force method used, unless
        gurobipy is available (default True)
    num_resamples: int
        The number of resamples used to compute the confidence intervals (default 1000)
    conf_level: float
        The confidence interval level (default 0.95)
    print_to_console: bool
        Print results directly to console (default False)
    seed: int
        Seed to generate a random number (default None)

    Returns
    -------
    time: numpy.array
        The Numpy array containing all the timepoints at which the simulation was run. The Numpy array is of dtype=float.
    mu: numpy.array
        The Numpy array containing the mean elementary effects over time. Each row is a timepoint and every column
        contains the mean elementary effects for a certain parameter. The Numpy array is of dtype=float.
    mu_star: numpy.array
        The Numpy array containing the absolute mean elementary effects over time. Each row is a timepoint and every column
        contains the absolute mean elementary effects for a certain parameter. The Numpy array is of dtype=float.
    sigma: numpy.array
        The Numpy array containing the standard deviation sof the mean elementary effect over time. Each row is a timepoint
        and every column contains the standard deviations of the mean elementary effect for a certain parameter.
        The Numpy array is of dtype=float.
    mu_star_conf_level: numpy.array
        The Numpy array containing the bootstrapped confidence intervals. Each row is a timepoint and every column
        contains the bootstrapped confidence interval for a certain parameter. The Numpy array is of dtype=float.
    """
    # First do some checks whether the problem variable contains all the needed keys for the function.
    if "num_vars" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'num_vars' (https://salib.readthedocs.io/en/latest/basics.html)")
    if "names" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'names' (https://salib.readthedocs.io/en/latest/basics.html)")
    if "bounds" not in problem:
        raise ValueError(
            "Variable 'problem' needs to contain key 'bounds' (https://salib.readthedocs.io/en/latest/basics.html)")

    # Determening the timepoints of the simulation.
    n = int(np.ceil(t_tot/dt) + 1)  # Number of timesteps of the simulation [-]
    time = np.linspace(0, t_tot, n)  # Array with all timepoints.

    # Retrieve number of parameters
    num_parameters = problem["num_vars"]

    # Generating input parameters for the model
    # Each column is a parameter, one row contains all input parameters for one simulation
    model_input = morris_sample.sample(
        problem, trajectories, num_levels=num_levels, optimal_trajectories=optimal_trajectories, local_optimization=local_optimization, seed=seed)

    # Check if you need to override the dna_conc variable then do np.delete()
    # https://stackoverflow.com/questions/24027040/how-to-extract-all-columns-but-one-from-an-array-or-matrix-in-python
    if "dna_conc" in problem["names"]:
        index = problem["names"].index("dna_conc")
        dna_conc_array = model_input[:, index]
        # Remove the column with the dna concentration, since this is a seperate value from the input parameters
        model_input_temp = np.delete(model_input, index, axis=1)
    else:
        dna_conc_array = np.ones(
            model_input.shape[0], dtype=np.float64) * dna_conc
        model_input_temp = model_input

    model_output = func(parameters=model_input_temp, dna_conc=dna_conc_array,
                        s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot)

    # Running the Morris analysis (using the output of all the different simulations)
    indices_dict = morris_analyze.analyze(
        problem, model_input, model_output, num_levels=num_levels)
    mu = indices_dict['mu']
    mu_star = indices_dict['mu_star']
    sigma = indices_dict['sigma']
    mu_star_conf_level = np.array(
        indices_dict['mu_star_conf'], dtype=np.float32)

    return mu, mu_star, sigma, mu_star_conf_level


def morris_problem_description_prokaryotic(s_i, problem, trajectories, num_levels, dt, t_tot, path="modelling/data", tag=f"_{int(time.time())}"):
    file_path = path + "/" + "description" + tag + ".txt"
    with open(file_path, "w") as wf:
        wf.write("Description Morris sensitivity analysis\n")
        wf.write("\n")
        wf.write(
            f"Morris sensitivity analysis for the prokaryotic system was run using {trajectories} trajectories and {num_levels} levels.\n")
        wf.write(
            f"The initial substrate concentration was {s_i} uM. The concentration of the product was taken as the output variable for the analysis.\n")
        wf.write(
            f"In total {len(problem['names'])} parameters were varied. The name of these parameters and their respective ranges are as follows:\n")
        wf.write("\n")
        count = 0
        for param in problem["names"]:
            # Write down the parameter name and the bounds
            wf.write(
                f"{param} = ({problem['bounds'][count][0]:.3e}, {problem['bounds'][count][1]:.3e})\n")
            # Next parameter is linked to next bound in the bounds list
            count += 1

        wf.write("\n")
        wf.write(
            f"The simulation was run for {t_tot} seconds with timesteps of {dt} seconds. Analysis was performed\n")
        wf.write("at every timestep.")


def morris_datawriter(problem, mu, mu_star, sigma, mu_star_conf_level, time=None, path=f"modelling/data", tag=f"_{int(time.time())}"):

    data_order = ["mu", "mu_star", "sigma", "mu_star_conf_level"]
    data_dict = {"mu": mu, "mu_star": mu_star, "sigma": sigma,
                 "mu_star_conf_level": mu_star_conf_level}

    # Only add time column if it is given
    if time is not None:
        fieldnames = ["time"]
    else:
        fieldnames = []

    # Create all the fieldnames
    for parameter in problem["names"]:
        fieldnames.append(parameter)

    for data_name in data_order:
        if "time" in fieldnames:
            data_to_write = np.concatenate(
                (time[:, np.newaxis], data_dict[data_name]), axis=1)
        else:
            data_to_write = data_dict[data_name]

        file_path = path + "/" + data_name + tag + ".csv"
        with open(file_path, "w") as wf:
            csv_writer = csv.writer(wf)
            csv_writer.writerow(fieldnames)

            if len(data_to_write.shape) == 1:
                csv_writer.writerow(data_to_write)
            elif len(data_to_write.shape) == 2:
                for i in range(data_to_write.shape[0]):
                    csv_writer.writerow(data_to_write[i, :])


def morris_datareader(path=f"modelling/data", tag=f"_{int(time.time())}"):

    data_names = ["mu", "mu_star", "sigma", "mu_star_conf_level"]
    data_dict = dict()

    for data_name in data_names:
        data = []
        file_path = path + "/" + data_name + tag + ".csv"
        with open(file_path, "r") as rf:
            csv_reader = csv.reader(rf)
            fieldnames = next(csv_reader)
            for row in csv_reader:
                data.append(row)

        data_array = np.array(data, dtype=np.float32)
        if "time" in fieldnames:
            parameters = fieldnames[1:]
            time_array = data_array[:, 0]
            data_array = np.delete(data_array, 0, axis=1)
        else:
            parameters = fieldnames
            time_array = None

        data_dict[data_name] = data_array

    # This only needs to be done once
    data_dict["time"] = time_array

    return parameters, data_dict


# def morris_datawriter(problem, path, filenumber, time, mu, mu_star, sigma, mu_star_conf_level):
#     # Names of the 4 files
#     filenames = ["\mu", "\mu_star", "\sigma", "\mu_star_conf_level"]
#     # Names of the columns of each file, namely time and all the parameters
#     fieldnames = ["time"]
#     for parameter in problem["names"]:
#         fieldnames.append(parameter)
#     # Put all data in 1 list
#     indices = [mu, mu_star, sigma, mu_star_conf_level]
#     # A column array containing all timepoints
#     time_column = np.zeros((time.shape[0], 1), dtype=np.float64)
#     time_column[:, 0] = time
#     # Loop over the four files
#     for ii in range(4):
#         # Make file
#         new_file_name = path + filenames[ii] + "_" + filenumber + ".csv"
#         new_file = open(new_file_name, "x")
#         csv_writer = csv.writer(new_file)
#         # Write header
#         csv_writer.writerow(fieldnames)
#         # Store all the information to be saved in the data variable
#         data = np.concatenate((time_column, indices[ii]), axis=1)
#         # Loop over rows of each file (each row denoting a timepoint)
#         for jj in range(time.shape[0]):
#             csv_writer.writerow(data[jj, :])
#         new_file.close()


# def morris_datareader(parameter, index, path, filenumber):
#     filename = path + "\\" + index + "_" + filenumber + ".csv"
#     print(filename)
#     file = open(filename, "r")
#     csv_reader = csv.DictReader(file)
#     data = []
#     for line in csv_reader:
#         data.append(line[parameter])
#     file.close()
#     map(int, data)
#     data = np.array(data, dtype=np.float64)
#     return data
