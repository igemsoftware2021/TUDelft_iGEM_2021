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
    sigma = data_dict["sigma"].reshape(-1)

    custom_colors = custom_aptavita_colors()
    factor_names = senstivity_analysis_factor_names()

    fig1, ax1 = plt.subplots(figsize=(14, 8), dpi=125)
    fig2, ax2 = plt.subplots(figsize=(14, 8), dpi=125)
    fig3, ax3 = plt.subplots(figsize=(14, 8), dpi=125)

    mu_fc = mu/standard_area
    mu_star_fc = mu_star/standard_area
    sigma_fc = sigma/standard_area

    ax1.bar(np.arange(0, mu_fc.shape[0]), mu_fc, color=custom_colors)
    ax2.bar(np.arange(0, mu_star_fc.shape[0]), mu_star_fc, color=custom_colors)
    ax3.bar(np.arange(0, sigma_fc.shape[0]), sigma_fc, color=custom_colors)

    # Set all proporties for ax1
    # ax1.set_xlabel(r"Time $[s]$")
    ax1.set_xticks(np.arange(0, len(parameters)))
    ax1.set_xticklabels(factor_names[:len(parameters)])
    ax1.set_ylabel(r"$\mathrm{{\mu}}$ $[\mathrm{\mu M \cdot s}]$")
    # ax1.yaxis.set_major_locator(MultipleLocator(5000))

    # Set all proporties for ax2
    ax2.set_xticks(np.arange(0, len(parameters)))
    ax2.set_xticklabels(factor_names[:len(parameters)])
    ax2.set_ylabel(r"$\mathrm{{\mu}}^{\ast}$ $[\mathrm{\mu M \cdot s}]$")
    # ax2.yaxis.set_major_locator(MultipleLocator(5000))

    # Set all proporites for ax3
    ax3.set_xticks(np.arange(0, len(parameters)))
    ax3.set_xticklabels(factor_names[:len(parameters)])
    ax3.set_ylabel(r"$\mathrm{{\sigma}}$ $[\mathrm{\mu M \cdot s}]$")
    # ax3.yaxis.set_major_locator(MultipleLocator(5000))

    if save_path is not None:
        fig1.savefig(f"{save_path}_mu_fc.svg", format="svg", dpi=1200)
        fig2.savefig(f"{save_path}_mu_star_fc.svg", format="svg", dpi=1200)
        fig3.savefig(f"{save_path}_sigma_fc.svg", format="svg", dpi=1200)
    else:
        plt.show()


def morris_method_visualization():
    num_factors = 3
    num_levels = 4
    num_trajectories = 5
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
    np.random.set_state(('MT19937', np.array([2147483648, 3518626697, 4058552728,  551327943, 4037528684,
                                              1617476776, 1582366576, 1127218356,  695230299,  887251314,
                                              1434336238,  562498140, 3919759145, 3667742909,  257426968,
                                              2558015267, 1393037283, 2782120227, 3761450755, 2047987210,
                                              3794524768, 4180794033,  866641479, 3679433940, 2177811301,
                                              554209624, 1673895573, 4144678097, 1929015793, 2820104439,
                                              789703441, 3893488458, 1806535622, 2282357178, 1177870788,
                                              1958531751, 3228327079, 3813284106, 3382292793, 2240830475,
                                              3414810018,  451715068, 1947396098, 4094714315,  757810080,
                                              807249853, 3602073203, 1416317232, 2880173888, 3416565716,
                                              3338684848,   85990715, 3647897893, 2537223738,  373833000,
                                              1692231494,  889132628,  974276034, 3173593593, 3210545286,
                                              3834922911, 1571823056, 1756660532, 3277489153, 2876771091,
                                              1838461604,  677162284, 3511393104, 1844015394,  488553694,
                                              1623690983,  956449945, 4098048616, 4061173888, 3910410884,
                                              768517141, 2180001659,  149137347,  845366307, 1982994281,
                                              2446249253, 1436529775, 3061520971, 2690839161, 2286555137,
                                              1957157927,  375911559, 1888506123, 4233428735, 1911289019,
                                              2872791403, 3220404776, 3680696425, 3819589194, 3703451070,
                                              146049519, 2910907883, 2287810285, 3338896371,  140987516,
                                              553540661, 3966161702, 3007575779, 1581569503, 4137590965,
                                              2489486962, 2125866918, 1896568528,  968963914, 2042184985,
                                              1239495713, 3120947550, 3160531082, 2959134992, 1000591379,
                                              1028950863,  740004944, 4183990114,  416981178,  421178497,
                                              3031717407, 3839635468,  599534188, 1343891905,  827036224,
                                              1216277790, 3024366788,  779676451, 1560894818, 1794244289,
                                              787213467,   46197467, 2573288874,  690311141,  793908338,
                                              3223822268,  236517427, 1537156111, 2518982324, 1836480173,
                                              1033401107,  950987918, 3874686387, 1055777206,  529810303,
                                              3340341737, 1570362537, 1166270872, 2039185856,  489631919,
                                              2620289998, 2568832646, 3772326514,  938438079,  913121597,
                                              764759630, 2944986243, 1692071482, 2594612390, 1105015430,
                                              538173560, 2333964461,  467073574, 3433454269, 2062427685,
                                              1316932158, 3543318930, 3990398019, 1773936535, 2409898164,
                                              701513664, 1956032845, 3182375716, 4127528869, 2014712534,
                                              3358360951,  381110876,  544306689,  302734011, 1906997681,
                                              2234359142,  941314843, 2305092612, 3654609487, 3314239264,
                                              20195556,  346264179, 3270710568, 3867485582, 2545941983,
                                              3113147406, 3048231659, 2338327307, 2060399615, 1270453449,
                                              249193665, 2962426968,  162079337, 1201467324, 1447463901,
                                              4016699714, 2692429035, 1491541367, 1542005298,  143648260,
                                              1845663522, 1176049241, 2305185853,  503128342, 1909701155,
                                              2646342539,   88049444, 1403939813, 3723718515, 1156274187,
                                              2652331916, 1330654715, 1593864272, 4018010925, 1238744078,
                                              1423906772,  519792923, 3888876747,  672192843,  247025362,
                                              1910462831, 1894714911, 2692846746, 1655275106, 2559791571,
                                              840945087, 3102425968, 2285977194, 2511655332,  708368680,
                                              473537506, 4163360428, 3450817295, 4201557422, 3444901390,
                                              3489721579,  816806879, 1582865384, 3347653010, 3428496401,
                                              4224592846, 1433380623,  648463165,  312692744, 1226323912,
                                              264244240, 2443458687, 3851224710, 1769711131, 2347276648,
                                              4144992852, 2338622223,   60423593, 1992985274, 3253963240,
                                              363039080,  631887209,  761452610, 3403564357,  170679303,
                                              228893748, 2466255248, 1535709266, 1622019621, 2005428648,
                                              752062867, 1736058913, 3583240568, 2788943563, 3417409665,
                                              24171883,  673951840,  298395209, 1120693036,   38895050,
                                              1959155309,  389157730, 2092199396, 1767914171,  300380908,
                                              1740192601,  804698102, 2622635785,  457030486, 1494025405,
                                              317008149,  545252964, 2138604467, 3128603706, 2723639031,
                                              1278227883,  221707037, 1671740084, 2401976123, 2317643775,
                                              3782094821,  129379350, 1941780349, 1364257803, 3366384717,
                                              2602228586, 2394729484,  385646311,  254719845,  708181813,
                                              3028405217,  967023160,  130614626, 3835508557, 1350955093,
                                              3618033789, 3227657065, 1739150281,  843550151,  500894336,
                                              3117996467, 2639343201, 2868097923, 3283573164, 3561785479,
                                              100380731,  454842684,  303652738, 1759772360, 2972585918,
                                              2543724558, 3446367408,  705045645, 1581128990, 4064393217,
                                              497380606,  860426146, 2833456403, 1685931085, 4215900498,
                                              1573561729, 1155266651, 2022788721, 1637257229, 4247461760,
                                              1936763424, 2021619247, 2554056160, 2434732564, 3328739208,
                                              2620873015, 4001287973,  115549832, 3115797255, 2354533543,
                                              543310969, 3435532301, 3462241565,  503678885, 3810506884,
                                              615256711, 2707607410,  409236852, 3406265564, 3862931636,
                                              846062285, 2974736526, 2805357677, 3938244183, 3713615913,
                                              1837240882, 2934548722, 3794598503, 3579383641, 3019423640,
                                              265378378,    4039318, 4156874516, 1012670860, 1237635103,
                                              3738033606, 1341380243,  943412411, 3409010282,  305918544,
                                              775284229, 2771293022,  938807278, 1455188599, 3248405342,
                                              2285635356, 2428858206, 2789281898,  135912329, 3127879423,
                                              211564359, 2786920598,  416379121, 4249176581, 3402458518,
                                              232001855, 1367137164,  806235454, 2412923849, 1349478436,
                                              2080324223, 1579818442, 2419933094, 2561677851, 2597336192,
                                              2666436030, 2717454882,  701016771, 1337223508, 3892057596,
                                              3960057682,  960928045, 4165525496, 2628708702,  615432885,
                                              855879624, 1225904465, 3202930760,  927674064, 1436287927,
                                              3942770705, 1951489601, 3690778151,    6745837,  159546402,
                                              2751830725,  857857658, 2420608172, 3399893206,  993978254,
                                              613268630, 4099640336,  856615950, 1510196067, 3787928156,
                                              4248018840,  125270695, 1913186520, 1158263506, 2614343896,
                                              3143255354, 1407404643, 1059305592, 1266887656, 3998192511,
                                              2118651192, 1158795018, 2463010309,  968754814, 3309312352,
                                              1235039238, 3990100991, 3880849345, 2470955735, 3355139616,
                                              275529447, 1870824098, 1238330767, 4203371473, 4086099815,
                                              1735563394, 2815712740, 2820016653, 1257165477, 2271166943,
                                              1894391818, 2495280417, 2608388563, 2671582317,  295188936,
                                              3773277500, 3328075093,  535893009, 2006730908, 3553890959,
                                              1833717031, 1253119489, 2657848475, 2229551359, 3988885770,
                                              680793556, 3488676944,  320264333,  176135252, 3630973052,
                                              3539023329, 4253964885, 2200702375, 1816688360, 3031670643,
                                              678118192,  418317316, 4263543120,  180433777, 2256384102,
                                              1302197718, 1330948948,  268777061, 1921107738, 1413183120,
                                              3172437959, 1381582217,  764339285, 3649763969, 3912469075,
                                              3678813638, 2595397609, 3617789896, 1297514342, 3858593165,
                                              512915081, 1162254733, 3174021615, 2969827363, 2498411377,
                                              2347417025, 1860866938, 2845064964, 2169313942, 2606978536,
                                              1401176484, 3222098253, 1164859957, 1283608681,  420900326,
                                              2305777247, 2550709267, 4257940864, 3204314493,  991905535,
                                              606878647, 1238125415, 3730466314,  447823841, 4151399594,
                                              1558749155, 3838448682, 2488006731, 1444275534, 1606475691,
                                              1286593823, 2689339324,  371720243,  628193722, 3497415222,
                                              3070793640,  575843273, 1648800877, 1840857257, 2481281013,
                                              170856436, 2220141279, 1740560414,  967057764, 1131466283,
                                              599998248,  835143095, 1808260471, 2707106112, 1335711702,
                                              1213702487, 3171708524, 2108296753, 1234275130, 2270691464,
                                              4074661390, 1168925303,  958691887, 2776395036, 3945429162,
                                              3314578123, 2840324384,  450756077, 2802322040,   73317957,
                                              640195299, 3306949852,    7509586, 1497085088, 3494544993,
                                              2010211649, 2826271717,  219477199, 1208690304, 3283547286,
                                              1405136256,  474942740, 3780779441, 1114024854, 3557711934,
                                              9814670, 2926482981, 1098879487, 2112419592,  352748701,
                                              1890898463,   55938565, 3285311873,  317548342,  201387787,
                                              3415547509, 3521136864,  579982282, 2411106269,  907735210,
                                              827555235, 1109587513, 1521239655, 2069216277, 1638501385,
                                              3526764973, 1348773455, 2984670198,  889304081, 3367575394,
                                              401914239, 2240995634, 2339147436, 2513608329], dtype=np.uint32), 623, 0, 0.0))

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
                        point[j-1] = point[j-1] - 1
                    else:
                        point[j-1] = point[j-1] + 1
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
    # morris_method_visualization()
    plot_morris_analysis_mu_star_subplots(
        path="modelling/data", tag="_1633520689", save_path="modelling/data/plots/T--TUDelft--Morris_Mu_Star_Subplots_1633520689.svg")
    # plot_morris_analysis_area(
    #     path="modelling/data", tag="_1633591400", save_path="modelling/data/plots/T--TUDelft--Morris_Area_1633591400")

    parameters = standard_parameters_prokaryotic()
    standard_area = model_prokaryotic_area(
        parameters, 3*10**-3, 250, 0.05, 0.09)
    plot_morris_analysis_area_fold_change(
        path="modelling/data", tag="_1633591400", save_path="modelling/data/plots/T--TUDelft--Morris_Area_1633591400", standard_area=standard_area)
