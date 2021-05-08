# Libraries to import

import numpy as np
import matplotlib.pyplot as plt
import math

# Constants

tau     = 0.63           # []
R       = 12*10**-6      # [m]
rho     = 1025           # [kg/m3]
K       = 5.92*10**-9    # [m3/m2/s]
µ       = 1.1*10**-3     # [kg/s^2/m]
theta   = math.pi/6      # [radians]
sigma   = 5.25*10**-2    # [kg/s^2]


# for loop for lenght L and length of tortuous hollow fiber lp

L = np.linspace(0, 10*10**-2, 101)    # from 0 to 10 cm, with steps of 1 mm [m] 

for x in L:
    lp  = tau * L                                                                     # [m]
    u   = ( 2 * L * sigma * math.cos(theta) * rho * K ) / ( lp * R * µ * (rho+L) )    # [m/s]

# Plotting of results

plt.plot(L, u)
plt.xlabel("Lenght L [m]")
plt.ylabel("Velocity u [m/s]")
plt.show()