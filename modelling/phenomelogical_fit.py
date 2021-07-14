import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

cwd = os.getcwd()
path = cwd + '/GFP_PURE.xlsx'

df = pd.read_excel (path, sheet_name=2)
data_GFP = np.array(df)


# Fitting function
def func(t, k_prime, k, K, n):
    return k_prime+k*((t**n)/(t**n+K**n))
# Experimental t and y data
x_Data = data_GFP[5, 3::] # Time in [h]
y_Data = data_GFP[8, 3::]

# Plot experimental data points
plt.plot(x_Data, y_Data, 'bo', label='experimental data')

# Initial guess for the parameters only if this is to be incorporated
# initialGuess = [1.0, 1.0]

# Perform the curve-fit
popt, pcov = curve_fit(func, x_Data, y_Data)
print(popt)


# Plot the fitted function
plt.plot(x_Data, func(x_Data, *popt), 'r', label='fitted parameters')