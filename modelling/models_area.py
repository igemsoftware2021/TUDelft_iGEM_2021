import numpy as np
from numba import njit, prange
from models import model_prokaryotic_readout


@njit(cache=True, nogil=True)
def model_prokaryotic_readout_area(parameters, constants, dna_conc, s_i, vit_conc1, vit_conc2, dt=0.01, t_tot=7200):
    initial_conditions1 = np.array([dna_conc, s_i, vit_conc1])
    initial_conditions2 = np.array([dna_conc, s_i, vit_conc2])

    time, absorbance1 = model_prokaryotic_readout(
        parameters, constants, initial_conditions1, dt=dt, t_tot=t_tot)
    _, absorbance2 = model_prokaryotic_readout(
        parameters, constants, initial_conditions2, dt=dt, t_tot=t_tot)

    # TODO determine whether area over time? (cumsum) or just total area?
    area = np.sum((absorbance2 - absorbance1)*dt)

    return area


@njit(parallel=True)
def model_prokaryotic_readout_area_parallel(parameters, constants, dna_conc, s_i, vit_conc1, vit_conc2, dt=0.01, t_tot=7200):
    # The amount of time steps
    n = int(np.ceil(t_tot/dt) + 1)
    # The number of simulations
    num_simulations = parameters.shape[0]

    # TODO if total area, then array should be smaller, just (n,1)
    model_output = np.zeros((n, num_simulations))
    # Every column is a unique simulation
    for ii in prange(num_simulations):
        _, model_output[:, ii] = model_prokaryotic_readout_area(
            parameters[ii, :], constants, dna_conc=dna_conc, s_i=s_i, vit_conc1=vit_conc1, vit_conc2=vit_conc2, dt=dt, t_tot=t_tot)
    return model_output
