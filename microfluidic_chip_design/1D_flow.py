# Libraries to import

import numpy as np
import matplotlib.pyplot as plt
import math

# Constants

tau     = 0.63           # []
R       = 12*10**-6      # [m]
K       = 6.39*10**-13   # [m^2]
µ       = 1.10*10**-3    # [kg/s^2/m]
theta   = math.pi/6      # [radians]
sigma   = 52.456*10**-3  # [N/m]
k       = 6.39*10**-13   # [m^2]


# for loop for lenght L and length of tortuous hollow fiber lp

L = np.linspace(0, 20*10**-2, 101)    # from 0 to 10 cm, with steps of 1 mm [cm] 

for x in L:
    lp  = (tau * L)                                                       # [cm]
    u   = (((2*sigma*math.cos(theta)*K)/(R*µ))*(k/((lp/L)*k+L*K)))*10**6  # [µm/s]

# Plotting of results

plt.plot(L, u)
plt.xlabel("Lenght L [mm]")
plt.ylabel("Velocity u [µm/s]")
plt.show()