import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

# Plot helpers


def custom_aptavita_colors():
    return ["#9B0138", "#057D54", "#FFCE3A", "#EDB4D7", "#4FD590", "#4D94EF", "#EFA54D", "#313175",
            "#667817", "#6FC0A8", "#D6682A", "#F3758A", "#755F26", "#A74D36", "#C88E99", "#A79536"]


def plot_ambient_light_different_conditions(save_path: str = None):

    # Data (structure [value sensor 1, value sensor 2, value sensor 3, value sensor 4])
    # Pure sensor values
    # mean_condition1 = np.array([393.4, 404.4, 398.9, 434.4])
    # mean_condition2 = np.array([82, 80.7, 91, 82.0])
    # mean_condition3 = np.array([6, 3, 4, 5])
    # mean_condition4 = np.array([0, 0, 0, 0])

    # sem_condition1 = np.array([0.0971, 0.163, 0.131, 0.223])
    # sem_condition2 = np.array([0, 0.0586, 0, 0.0277])
    # sem_condition3 = np.array([0, 0, 0, 0])
    # sem_condition4 = np.array([0, 0, 0, 0])

    # Transmitted light normalization
    mean_condition1 = np.array([100, 100, 100, 100])
    mean_condition2 = np.array([20.84, 19.96, 22.82, 18.87])
    mean_condition3 = np.array([7.317, 3.717, 4.396, 6.101])
    mean_condition4 = np.array([0, 0, 0, 0])

    conf_sem_condition1 = np.array([0.0967, 0.119, 0.0970, 0.152])
    conf_sem_condition2 = np.array([0.01, 0.0442, 0.0147, 0.0315])
    conf_sem_condition3 = np.array([0.00354, 0.00293, 0.00283, 0.00614])
    conf_sem_condition4 = np.array([0, 0, 0, 0])

    custom_colors = custom_aptavita_colors()

    fig1, ax1 = plt.subplots(figsize=(12, 6), dpi=150)

    # Position of bars on x-axis
    ind = np.arange(mean_condition1.shape[0])

    # Width of a bar
    width = 0.20

    ax1.bar(ind-1.5*width, mean_condition1, width,
            color=custom_colors[0], label="Condition 1")
    ax1.bar(ind-0.5*width, mean_condition2, width,
            color=custom_colors[1], label="Condition 2")
    ax1.bar(ind+0.5*width, mean_condition3, width,
            color=custom_colors[2], label="Condition 3")
    ax1.bar(ind+1.5*width, mean_condition4, width,
            color=custom_colors[3], label="Condition 4")

    # Set all proporties for ax1
    ax1.legend()
    ax1.set_xticks(np.arange(0, 4))
    ax1.set_xticklabels(["Sensor 1", "Sensor 2", "Sensor 3", "Sensor 4"])
    # ax1.set_yscale("log")
    ax1.set_ylabel("Transmitted light (%)")

    ax1.yaxis.set_major_locator(MultipleLocator(10))
    ax1.yaxis.set_minor_locator(MultipleLocator(5))

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def absorbance_spectrum_food_colorant(save_path: str = None):

    wavelength = np.array([300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630,
                          640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000])
    absorbance = np.array([1.237300038, 1.203600049, 1.036299944, 0.6395999789, 0.3312999904, 0.244599998, 0.2535000145, 0.3095000088, 0.4287, 0.5939000249, 0.7129999995, 0.9079999924, 0.8080000281, 0.4555000067, 0.181400001, 0.0903000012, 0.07079999894, 0.07079999894, 0.07859999686, 0.09089999646, 0.1116999984, 0.1460999995, 0.1965000033, 0.2802999914, 0.4092999995, 0.5819000006, 0.8377000093, 1.20630002, 1.74000001, 2.361299992, 2.769500017, 3.173199892, 3.296099901, 3.650399923, 3.471400023, 3.382100105,
                          2.788899899, 1.676300049, 0.8014000058, 0.3761999905, 0.1768999994, 0.09860000014, 0.06849999726, 0.05469999835, 0.04850000143, 0.04600000009, 0.04450000077, 0.04369999841, 0.04300000146, 0.0428000018, 0.04239999875, 0.04230000079, 0.04250000045, 0.04230000079, 0.04210000113, 0.04210000113, 0.04259999841, 0.0447999984, 0.0449000001, 0.04300000146, 0.04259999841, 0.04300000146, 0.04390000179, 0.04410000145, 0.04390000179, 0.0438000001, 0.04399999976, 0.04450000077, 0.04580000043, 0.04569999874, 0.04580000043])
    absorbance_norm = absorbance / np.amax(absorbance)

    fig1, ax1 = plt.subplots(figsize=(10, 5), dpi=150)
    ax1.plot(wavelength, absorbance_norm, color="#9B0138")
    ax1.set_xlabel("Wavelength [nM]")
    ax1.set_ylabel("Absorbance")

    ax1_ylim = ax1.get_ylim()
    ax1.set_ylim((0, ax1_ylim[1]))

    ax1.xaxis.set_major_locator(MultipleLocator(50))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.yaxis.set_major_locator(MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.01))

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_absorbance_spectrum_cpr(save_path: str = None):

    wavelength = np.array([300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630,
                          640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000])
    absorbance = np.array([2.1627, 1.9501, 1.8003, 1.716, 1.6937, 1.6973, 1.715, 1.7267, 1.7012, 1.6293, 1.5652, 1.5446, 1.5455, 1.5574, 1.5718, 1.5872, 1.6043, 1.6281, 1.66, 1.6962, 1.7507, 1.8194, 1.9025, 2.0061, 2.0973, 2.194, 2.3261, 2.4813, 2.4542, 2.2637, 1.9338, 1.6342, 1.5003, 1.4298,
                          1.4075, 1.3963, 1.3881, 1.3835, 1.3828, 1.3807, 1.3774, 1.3765, 1.3734, 1.3722, 1.3697, 1.368, 1.3677, 1.3645, 1.3619, 1.3616, 1.3598, 1.3575, 1.3566, 1.3561, 1.354, 1.352, 1.3515, 1.3519, 1.3474, 1.3466, 1.3458, 1.3431, 1.3435, 1.3402, 1.3377, 1.3374, 1.3362, 1.3351, 1.3344, 1.3307, 1.3303])
    absorbance_blanco = np.array([1.7635, 1.623, 1.5849, 1.5475, 1.5199, 1.4977, 1.48, 1.4671, 1.4518, 1.4439, 1.4358, 1.4353, 1.4316, 1.4279, 1.4273, 1.424, 1.4218, 1.4197, 1.4158, 1.4161, 1.4143, 1.4123, 1.4096, 1.4064, 1.4058, 1.4047, 1.4035, 1.4026, 1.404, 1.399, 1.3952, 1.3942, 1.3919, 1.3926,
                                 1.389, 1.3879, 1.3841, 1.3834, 1.3828, 1.3785, 1.3761, 1.3755, 1.3737, 1.3711, 1.3689, 1.3685, 1.366, 1.3638, 1.3606, 1.361, 1.3598, 1.3582, 1.3547, 1.3569, 1.3547, 1.3524, 1.3499, 1.352, 1.3498, 1.3473, 1.3444, 1.3438, 1.3445, 1.3399, 1.3386, 1.3385, 1.3377, 1.3362, 1.3364, 1.3316, 1.3299])

    absorbance = absorbance - absorbance_blanco

    absorbance_norm = absorbance / np.amax(absorbance)

    fig1, ax1 = plt.subplots(figsize=(12, 6), dpi=150)
    ax1.plot(wavelength, absorbance_norm, color="#9B0138")
    ax1.set_xlabel("Wavelength [nM]")
    ax1.set_ylabel("Absorbance")

    print(
        f"Wavelength maximum absorbance: {wavelength[np.argmax(absorbance_norm)]}")

    # ax1_ylim = ax1.get_ylim()
    # ax1.set_ylim((0, ax1_ylim[1]))

    ax1.xaxis.set_major_locator(MultipleLocator(50))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))

    ax1.yaxis.set_major_locator(MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.01))

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


def plot_absorbance_curves_hardware_spectr(save_path: str = None):

    cpr_conc_hardware = np.array([0.25, 0.25, 0.25, 0.5, 0.5, 0.5,
                                  0.75, 0.75, 0.75, 1, 1, 1, 1.25, 1.25, 1.25])
    cpr_absorbance_hardware = np.array([0.29901, 0.24646, 0.27459, 0.62725, 0.577777, 0.67268067, 0.8091459787, 0.7160924046,
                                        0.6714569699, 0.9095247866, 0.9010765457, 0.9437113038, 1.0697792, 1.115627537, 1.18410570])

    cpr_conc_spectr = np.array(
        [0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1.25, 1.25])
    cpr_absorbance_spectr = np.array(
        [0.3045, 0.2613, 0.53305, 0.50545, 0.60035, 0.62765, 0.78615, 0.70985, 0.95405, 0.90065])

    # First all calculations
    def func(x, a, b):
        return a*x + b

    # wh = without
    def func_wh_offset(x, a):
        return a*x

    # Hardware calculations
    popt_func_hardware, pcov_func_hardware = curve_fit(
        func, cpr_conc_hardware, cpr_absorbance_hardware)
    popt_func_wh_offset_hardware, pcov_func_wh_offset_hardware = curve_fit(
        func_wh_offset, cpr_conc_hardware, cpr_absorbance_hardware)

    cpr_conc_hardware_line = np.arange(0.05,
                                       np.amax(cpr_conc_hardware)+0.1, step=0.05)

    cpr_absorbance_func_hardware = func(
        cpr_conc_hardware_line, popt_func_hardware[0], popt_func_hardware[1])

    cpr_absorbance_func_wh_offset_hardware = func_wh_offset(
        cpr_conc_hardware_line, popt_func_wh_offset_hardware[0])

    # Calculate R^2 func hardware
    residuals_func_hardware = cpr_absorbance_hardware - \
        func(cpr_conc_hardware, popt_func_hardware[0], popt_func_hardware[1])
    ss_res_func_hardware = np.sum(residuals_func_hardware**2)
    ss_tot_func_hardware = np.sum(
        (cpr_absorbance_hardware-np.mean(cpr_absorbance_hardware))**2)
    r_squared_func_hardware = 1 - (ss_res_func_hardware/ss_tot_func_hardware)

    # Calculate R^2 func_wh_offset hardware
    residuals_func_wh_offset_hardware = cpr_absorbance_hardware - \
        func_wh_offset(cpr_conc_hardware, popt_func_hardware[0])
    ss_res_func_wh_offset_hardware = np.sum(
        residuals_func_wh_offset_hardware**2)
    ss_tot_func_wh_offset_hardware = np.sum(
        (cpr_absorbance_hardware-np.mean(cpr_absorbance_hardware))**2)
    r_squared_func_wh_offset_hardware = 1 - \
        (ss_res_func_wh_offset_hardware/ss_tot_func_wh_offset_hardware)
    #
    #
    #
    #
    # Spectrophotometer
    popt_func_spectr, pcov_func_spectr = curve_fit(
        func, cpr_conc_spectr, cpr_absorbance_spectr)
    popt_func_wh_offset_spectr, pcov_func_wh_offset_spectr = curve_fit(
        func_wh_offset, cpr_conc_spectr, cpr_absorbance_spectr)
    print(f"popt hardware: ", popt_func_hardware)
    print(f"popt spectr: ", popt_func_spectr)

    cpr_conc_spectr_line = np.arange(
        0.05, np.amax(cpr_conc_spectr)+0.1, step=0.05)

    cpr_absorbance_func_spectr = func(
        cpr_conc_spectr_line, popt_func_spectr[0], popt_func_spectr[1])

    cpr_absorbance_func_wh_offset_spectr = func_wh_offset(
        cpr_conc_spectr_line, popt_func_wh_offset_spectr[0])

    # Calculate R^2 func spectrophotometer
    residuals_func_spectr = cpr_absorbance_spectr - \
        func(cpr_conc_spectr, popt_func_spectr[0], popt_func_spectr[1])
    ss_res_func_spectr = np.sum(residuals_func_spectr**2)
    ss_tot_func_spectr = np.sum(
        (cpr_absorbance_spectr-np.mean(cpr_absorbance_spectr))**2)
    r_squared_func_spectr = 1 - (ss_res_func_spectr/ss_tot_func_spectr)

    # Calculate R^2 func_wh_offset spectrophotometer
    residuals_func_wh_offset_spectr = cpr_absorbance_spectr - \
        func_wh_offset(cpr_conc_spectr, popt_func_spectr[0])
    ss_res_func_wh_offset_spectr = np.sum(
        residuals_func_wh_offset_spectr**2)
    ss_tot_func_wh_offset_spectr = np.sum(
        (cpr_absorbance_spectr-np.mean(cpr_absorbance_spectr))**2)
    r_squared_func_wh_offset_spectr = 1 - \
        (ss_res_func_wh_offset_spectr/ss_tot_func_wh_offset_spectr)

    print(
        f"R^2 func hardware: {r_squared_func_hardware}, R^2 func_wh_offset hardware: {r_squared_func_wh_offset_hardware}")
    print(
        f"R^2 func spectr: {r_squared_func_spectr}, R^2 func_wh_offset spectr: {r_squared_func_wh_offset_spectr}")

    fig1, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5), dpi=150)
    # Hardware plotting
    ax1.scatter(cpr_conc_hardware, cpr_absorbance_hardware,
                color="#057D54", alpha=1.0, s=15)
    ax1.plot(cpr_conc_hardware_line, cpr_absorbance_func_hardware,
             color="#9B0138", label="best fit")
    # ax1.plot(cpr_conc_hardware_line,
    #          cpr_absorbance_func_wh_offset_hardware, color="#FFCE3A", label="0 constrain")

    # Spectrophotometer plotting
    ax2.scatter(cpr_conc_spectr, cpr_absorbance_spectr,
                color="#057D54", alpha=1.0, s=15)
    ax2.plot(cpr_conc_spectr_line, cpr_absorbance_func_spectr,
             color="#9B0138", label="best fit")  # "#000000"
    # ax2.plot(cpr_conc_spectr_line,
    #          cpr_absorbance_func_wh_offset_spectr, color="#FFCE3A", label="0 constrain")

    # Proporties ax1
    ax1.set_xlabel("CPR concentration [mM]")
    ax1.set_ylabel("Absorbance")

    ax1_xlim = ax1.get_xlim()
    ax1.set_xlim((0, ax1_xlim[1]))
    ax1_ylim = ax1.get_ylim()
    ax1.set_ylim((0, ax1_ylim[1]))

    ax1.xaxis.set_major_locator(MultipleLocator(0.2))
    ax1.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax1.yaxis.set_major_locator(MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.05))

    # Proporties ax2
    ax2.set_xlabel("CPR concentration $[\mathrm{mM}]$")
    ax2.set_ylabel("Absorbance")

    ax2_xlim = ax2.get_xlim()
    ax2.set_xlim((0, ax2_xlim[1]))
    ax2_ylim = ax1.get_ylim()
    ax2.set_ylim((0, ax2_ylim[1]))

    ax2.xaxis.set_major_locator(MultipleLocator(0.2))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax2.yaxis.set_major_locator(MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.05))

    # Set the character labels
    ax1.text(-0.1, 1.05, "a", transform=ax1.transAxes,
             size=16, weight="bold")
    ax2.text(-0.1, 1.05, "b", transform=ax2.transAxes,
             size=16, weight="bold")

    if save_path is not None:
        fig1.savefig(f"{save_path}", format="svg", dpi=1200)
    else:
        plt.show()


if __name__ == "__main__":
    pass
    plot_ambient_light_different_conditions(
        "T--TUDelft--Hardware_Different_Light_Conditions.svg")
    # absorbance_spectrum_food_colorant()
    plot_absorbance_spectrum_cpr("T--TUDelft--CPR_Absorbance_Spectrum.svg")
    plot_absorbance_curves_hardware_spectr(
        "T--TUDelft--CPR_Absorbance_Curves_Hardware_Spectr.svg")
