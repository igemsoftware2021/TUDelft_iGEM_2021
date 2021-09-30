import numpy as np
import matplotlib.pyplot as plt

from morris_method import morris_datareader

parameters, data_dict = morris_datareader(
    path="modelling/data", tag="_1633001637")

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()

for i in range(len(parameters)):
    ax1.plot(data_dict["time"], data_dict["mu"][:, i], label=parameters[i])
    ax2.plot(data_dict["time"], data_dict["mu_star"]
             [:, i], label=parameters[i])
    ax3.plot(data_dict["time"], data_dict["sigma"][:, i], label=parameters[i])
    ax4.plot(data_dict["time"], data_dict["mu_star_conf_level"]
             [:, i], label=parameters[i])

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()


plt.show()
