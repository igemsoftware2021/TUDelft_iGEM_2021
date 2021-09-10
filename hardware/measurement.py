import time
import multiprocessing
import numpy as np
import pigpio
from helpers import initialize_light_pins, initialize_heating_pin, main_measurement
from helpers import i2c_multiplexer_select_channel, i2c_activate_als_all_sensors
from helpers import i2c_write_to_all_sensors, pre_heater
from helpers import save_as_csv, read_from_csv
from apds9930 import APDS9930
from apds9930.values import APDS9930_ATIME
import matplotlib.pyplot as plt

###############################################################################
# Set-up

# Constants for test
total_time = 4

# Control pins
pin_heating = 18
pin_light_1 = 17
# pin_light_2 = "TBD"
# pin_light_3 = "TBD"
# pin_light_4 = "TBD"
pins_light = [pin_light_1]  # , pin_light_2, pin_light_3, pin_light_4]

# Constants for heating
k_u = 1  # Ultimate gain, used for Ziegler-Nichols tuning of the PID controller
t_u = 1  # Oscillation frequency associated with the ultimate gain
k_p = 0.20 * k_u        # Proportional gain PID controller
k_i = 0.50 * t_u        # Integral gain PID controller
k_d = 0.40 * k_u / t_u  # Derivative gain PID controller
temperature_desired = 25  # Desired temperature in degrees Celsius
pwm_freq = 10  # simulations suggest this is fine, see LT spice. Need to do this properly with a good toggle speed
pid_parameters = [k_p, k_i, k_d, temperature_desired, pwm_freq]
v_ref = 3.3
gain = 2.96


# SPI
spi_channel = 0  # SPI channel of the Raspberry Pi that is connected to the ADC
spi_baud = 50000  # Speed of the serial communication in bits/s
spi_flags = 0

# I2C
i2c_bus = 1  # I2C bus of the Raspberry Pi that is connected to the multiplexer
i2c_multiplexer_adress = 0x70  # TCA9548
i2c_sensor_adress = 0x39  # APDS9930
i2c_sensor_int_time = 0xc0  # Integration time APDS9930. 0x0c <-> 175 ms
i2c_sensor_gain = "TBD"  # Gain of the APDS9930.
i2c_sensor_1_channel = 0  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_2_channel = 1  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_3_channel = 2  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_4_channel = 3  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_channels = [i2c_sensor_1_channel]  # , i2c_sensor_2_channel,
# i2c_sensor_3_channel, i2c_sensor_4_channel]
###############################################################################


# Initialization

# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()
    # TODO display error

# Initialize pins
initialize_light_pins(pi, pins_light)
initialize_heating_pin(pi, pin_heating)


# Initialize SPI connection
adc_handle = pi.spi_open(spi_channel, spi_baud, spi_flags)

# Initialize I2C connection
i2c_multiplexer_handle = pi.i2c_open(i2c_bus, i2c_multiplexer_adress)
i2c_multiplexer_select_channel(pi,
                               i2c_multiplexer_handle, i2c_sensor_1_channel)  # Switch to a sensor channel
i2c_sensor_handle = APDS9930(i2c_bus)  # Create 1 sensor handle

# Activate ALS, and set gain and integration time
i2c_activate_als_all_sensors(
    pi, i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels)
i2c_write_to_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels,
                         APDS9930_ATIME, i2c_sensor_int_time)
# helpers.i2c_change_gain_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, gain)

# TODO test gain

###############################################################################


# print("Would you like to pre-heat the box? (yes or no)")
# pre_heating = input()
# if pre_heating == "yes":
#     pre_heater(pi, adc_handle, spi_channel, pin_heating,
#                pid_parameters, v_ref, gain)
#     print("Pre-heating has finished.")

# Heat

# print("Would you like to run a test? (yes or no)")
# test = input()
test = "yes"
if test == "yes":
    path = "hardware/test_results/"
    name = ["sensor_1", "sensor_2", "sensor_3", "sensor_4"]

    timepoints, absorbance, temperature_error = main_measurement(pi, pins_light, pin_heating,
                                                                 i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels, adc_handle,
                                                                 spi_channel, total_time, pid_parameters, v_ref, gain,
                                                                 duration=0.10, interval=0.005)  # returns a list of np arrays
    print(temperature_error)
    print(timepoints[-1])
    for i in range(len(i2c_sensor_channels)):
        timepoints_single_channel = timepoints[i]
        absorbance_single_channel = absorbance[i]
        save_as_csv(timepoints_single_channel,
                    absorbance_single_channel, path, name[i])
    save_as_csv(timepoints[-1], temperature_error, path, "temperature_error")
print("The test has been completed.")
timepoints_plot, absorbance_plot = read_from_csv(path, "sensor_1")
fig1, ax1 = plt.subplots()
ax1.plot(timepoints_plot, absorbance_plot)
fig1.show()

timepoints_plot, temperature_error = read_from_csv(path, "temperature_error")
fig2, ax2 = plt.subplots()
ax2.plot(timepoints_plot, temperature_error)
fig2.show()

plt.show()
