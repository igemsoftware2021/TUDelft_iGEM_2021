import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from tqdm import tqdm
from models import model_prokaryotic_hp
from standard_values import standard_parameters_prokaryotic, standard_constants, standard_initial_conditions


def plot_different_dt_error(low=0, high=-8, num=17, base=10, t_tot=7200):
    all_dt = np.logspace(low, high, num=num, base=base)

    parameters = standard_parameters_prokaryotic()
    constants = standard_constants()
    initial_conditions = standard_initial_conditions()

    time_final = np.arange(0, t_tot, 0.1, dtype=np.float32)
    p_final = np.zeros((num, time_final.shape[0]), dtype=np.float32)
    p_error = np.zeros((num, time_final.shape[0]), dtype=np.float32)

    for i in tqdm(range(all_dt.shape[0])):
        time, p = model_prokaryotic_hp(
            parameters, constants, initial_conditions, dt=all_dt[i], t_tot=t_tot, dt_store=0.01)

        f = UnivariateSpline(time, p, k=1, s=0, ext=2)
        p_final[i, :] = f(time_final)

    p_error = p_final - p_final[-1, :]
    p_error = p_error / p_final[-1, :]

    # Shorten the arrays due to some weird artifacts at the end
    time_final = time_final[:int(np.ceil(p_error.shape[1]*0.98))]
    p_error = p_error[:, :int(np.ceil(p_error.shape[1]*0.98))]

    fig, ax = plt.subplots()

    for i in range(p_error.shape[0]):
        ax.plot(time_final, p_error[i, :], label=f"{all_dt[i]}")
    ax.legend()
    # ax.set_yscale("log")
    # ax.set_xlim(0, t_tot*0.98)
    plt.show()


if __name__ == "__main__":
    plot_different_dt_error(low=0, high=-6, num=7, base=10, t_tot=7200)
