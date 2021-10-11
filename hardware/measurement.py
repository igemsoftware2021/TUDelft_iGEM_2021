import time
import multiprocessing
import numpy as np
import pigpio
from helpers import initialize_light_pins, initialize_heating_pin
from helpers import i2c_multiplexer_select_channel, i2c_activate_als_all_sensors
from helpers import i2c_write_to_all_sensors, pre_heater, temperature_controller
from helpers import determine_intensity_single_channel, calculate_absorbance
from helpers import save_as_csv, read_from_csv
from apds9930 import APDS9930
from apds9930.values import APDS9930_ATIME
import matplotlib.pyplot as plt
from file_helpers import write_temperature_csv, write_absorbance_csv
from file_helpers import read_absorbance_csv, read_temperature_csv
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
duration = 0.10
interval = 0.005

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
num_sensors = len(i2c_sensor_channels)
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


print("Pre-heat the box? (y or n)")
pre_heating = input()
if pre_heating == "y":
    temperature_time, temperature_error = pre_heater(pi, adc_handle, spi_channel, pin_heating,
                                                     pid_parameters, v_ref, gain)
    write_temperature_csv(temperature_time, temperature_error)


print("start test? (y or n)")
start_test = input()
if start_test == "y":
    start_time = time.time()
    stop_time = start_time + total_time
    intensity = []
    absorbance = []
    timepoints = []
    temperature_error = []

    # Add a list in every list, where every list is a channel
    for i in range(num_sensors):
        intensity.append([])
        absorbance.append([])
        timepoints.append([])
    timepoints.append([])

    T_c = 0.001  # Lowest high or low time [ms] to protect the hardware
    duty_cycle_lower_bound = T_c * pwm_freq * 10**6
    duty_cycle_upper_bound = (1 - T_c * pwm_freq) * 10**6

    print("start")
    while time.time() < stop_time:
        for i in range(num_sensors):
            pi.write(pins_light[i], 1)
            temperature_error = temperature_controller(pi, adc_handle, spi_channel, pin_heating, duty_cycle_lower_bound, duty_cycle_upper_bound,
                                                       pid_parameters, v_ref, gain, temperature_error, duration=duration, interval=interval)
            temperature_timepoint = time.time()
            timepoint, intensity_datapoint = determine_intensity_single_channel(pi,
                                                                                pins_light[i], i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers[i])
            intensity[i].append(intensity_datapoint)
            timepoints[i].append(timepoint - start_time)
            timepoints[-1].append(temperature_timepoint - start_time)

    pi.hardware_PWM(pin_heating, pwm_freq, 0)  # Turn of the PWM output signal
    print("test completed, processing results")
    for i in range(num_sensors):  # check this
        intensity_single_channel = intensity[i]
        absorbance[i] = calculate_absorbance(intensity_single_channel)

write_absorbance_csv(timepoints[0:num_sensors], absorbance)
write_temperature_csv(timepoints[-1], temperature_error)

print("finished processing result")
