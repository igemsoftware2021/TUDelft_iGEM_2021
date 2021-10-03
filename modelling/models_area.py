import numpy as np
from numba import njit, prange
from models import model_prokaryotic, model_prokaryotic_absorbance
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions
import matplotlib.pyplot as plt


@njit(cache=True, nogil=True)
def model_prokaryotic_area(parameters, dna_conc, s_i, vit_conc1, vit_conc2, dt=0.01, t_tot=7200):
    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])

    time, p1 = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions1, dt=dt, t_tot=t_tot)
    time, p2 = model_prokaryotic(
        parameters=parameters, initial_conditions=initial_conditions2, dt=dt, t_tot=t_tot)

    area = np.sum((p1 - p2)*dt)

    return area


@njit(parallel=True)
def model_prokaryotic_area_parallel(parameters, dna_conc, s_i, vit_conc1, vit_conc2, dt=0.01, t_tot=7200):
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]

    # TODO if total area, then array should be smaller, just (n,1)
    model_output = np.zeros(num_simulations, dtype=np.float64)
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        model_output[ii] = model_prokaryotic_area(
            parameters=parameters[ii, :], dna_conc=dna_conc[ii], s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot)
    return model_output


@njit(cache=True, nogil=True)
def model_prokaryotic_absorbance_area(parameters, constants, dna_conc, s_i, vit_conc1, vit_conc2, dt=0.01, t_tot=7200):
    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])

    time, absorbance1 = model_prokaryotic_absorbance(
        parameters=parameters, initial_conditions=initial_conditions1, constants=constants, dt=dt, t_tot=t_tot)
    time, absorbance2 = model_prokaryotic_absorbance(
        parameters=parameters, initial_conditions=initial_conditions2, constants=constants, dt=dt, t_tot=t_tot)

    area = np.sum((absorbance1 - absorbance2)*dt)

    return area


@njit(parallel=True)
def model_prokaryotic_absorbance_area_parallel(parameters, constants, dna_conc, s_i, vit_conc1, vit_conc2, dt=0.01, t_tot=7200):
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]

    # TODO if total area, then array should be smaller, just (n,1)
    model_output = np.zeros(num_simulations, dtype=np.float64)
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        model_output[ii] = model_prokaryotic_absorbance_area(
            parameters=parameters[ii, :], constants=constants, dna_conc=dna_conc[ii], s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot)
    return model_output


if __name__ == "__main__":
    parameters = standard_parameters_prokaryotic()
    # parameters2 = standard_parameters_prokaryotic()
    # parameters2[13] = parameters2[13]/10
    # parameters = np.vstack((parameters, parameters2))
    # print(parameters, parameters.shape)
    constants = standard_constants()
    # result = model_prokaryotic_absorbance_area(
    #     parameters, constants, 2*10**-3, 250, 0.05, 0.09)
    result = model_prokaryotic_area(
        parameters, 3*10**-3, 250, 0.05, 0.09)
    print(result)
    vit_conc = np.linspace(1, 20, 21)
    results = np.zeros(vit_conc.shape)
    # for i in range(vit_conc.shape[0]):
    #     results[i] = model_prokaryotic_readout_area(
    #         parameters, constants, 2*10**-3, 150, vit_conc[0], vit_conc[i], dt=0.01)

    # plt.plot(vit_conc, results)
    # plt.show()
