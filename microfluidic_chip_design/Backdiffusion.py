# Libraries to import

import numpy as np
import matplotlib.pyplot as plt
import math

# Constants

kB  = 1.38064852*10**-23 # [m2kg/s2K]
T   = 310.15             # [K]
pi  = math.pi            # [-]
eta = 1.10*10**-3        # [kg/ms]
r1  = 6.2*10**-9         # Radius for beta-galactosidase [m]
r2  = 4.2*10**-9         # Radius for catechol 2,3 dioxygenase [m]
r3  = 1.72*10**-9           # Estimated radius for CRDG [m]
r4  = 1.9*10**-9         # Estimated radius for pyrocathecol [m]


# Formula for diffusion coefficient D both enzymes

D1 = (kB*T)/(6*pi*eta*r1)  # D for beta-galactosidase [m2/s]

D2 = (kB*T)/(6*pi*eta*r2)  # D for catechol 2,3 dioxygenase [m2/s]



# for loop for time t and distance x both enzymes

t = np.linspace(0, 3600, 61)    # from 0 to 1 hour, with steps of 1 min [s] 

for i in t:
    x1 = np.sqrt(2*D1*t)  # Diffused distance for beta-galactosidase [s]
    x2 = np.sqrt(2*D2*t)  # Diffused distance for catechol 2,3 dioxygenase [s]

# Plotting of results both enzymes

plt.plot(t, x1, label="Beta-galactosidase")
plt.plot(t, x2, label="Catechol 2,3 dioxygenase")
plt.xlabel("Time t [s]")
plt.ylabel("Distance x [m]")
plt.legend()
plt.show()



# Formula for diffusion coefficient D both substrates

D3 = (kB*T)/(6*pi*eta*r3)  # D for CRDG

D4 = (kB*T)/(6*pi*eta*r4)  # D for pyrocatechol

# for loop for time t and distance x for both substrates

for i in t:
    x3 = np.sqrt(2*D3*t)  # Diffused distance for CRDG
    x4 = np.sqrt(2*D4*t)  # Diffused distance for pyrocatechol
    
# Plotting of results both substrates

plt.plot(t, x3, color='red', label="CRDG")
plt.plot(t, x4, color='black', label="Pyrocatechol")
plt.xlabel("Time t [s]")
plt.ylabel("Distance x [m]")
plt.legend()
plt.show()