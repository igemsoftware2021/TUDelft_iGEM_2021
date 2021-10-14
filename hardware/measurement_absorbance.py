import time
import multiprocessing
import numpy as np
import pigpio
from helpers import initialize_light_pins, initialize_heating_pin
from helpers import i2c_multiplexer_select_channel, i2c_activate_als_all_sensors
from helpers import i2c_write_to_all_sensors, pre_heater, temperature_controller
from helpers import determine_intensity_single_channel, calculate_absorbance
from helpers import i2c_set_gain_all_sensors
from apds9930 import APDS9930
from apds9930.values import APDS9930_ATIME
import matplotlib.pyplot as plt
from file_helpers import write_temperature_csv, write_absorbance_csv
from file_helpers import read_absorbance_csv, read_temperature_csv
###############################################################################
# Set-up

# Constants for test
total_time = 20

# Control pins
pin_light_1 = 17
pin_light_2 = 27
pin_light_3 = 22
pin_light_4 = 23
pins_light = [pin_light_1, pin_light_2, pin_light_3, pin_light_4]
# I2C
i2c_bus = 1  # I2C bus of the Raspberry Pi that is connected to the multiplexer
i2c_multiplexer_adress = 0x70  # TCA9548
i2c_sensor_adress = 0x39  # APDS9930
i2c_sensor_int_time = 0xc0  # Integration time APDS9930. 0x0c <-> 175 ms
i2c_sensor_gain = 2  # Gain of the APDS9930.
i2c_sensor_1_channel = 4  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_2_channel = 5  # Channel of the multiplexer to which the sensor is connected
# Channel of the multiplexer to which the sensor is connected
i2c_sensor_3_channel = 6
i2c_sensor_4_channel = 7  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_channels = [i2c_sensor_1_channel, i2c_sensor_2_channel,
                       i2c_sensor_3_channel, i2c_sensor_4_channel]
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
#helpers.i2c_change_gain_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, gain)
#i2c_sensor_handle.ambient_light_gain = 1

# TODO test gain
i2c_set_gain_all_sensors(pi, i2c_multiplexer_handle,
                         i2c_sensor_handle, i2c_sensor_channels, i2c_sensor_gain)
###############################################################################


start_test = "y"
if start_test == "y":
    start_time = time.time()
    stop_time = start_time + total_time
    intensity = []
    absorbance = []
    timepoints = []

    # Add a list in every list, where every list is a channel
    for i in range(num_sensors):
        intensity.append([])
        absorbance.append([])
        timepoints.append([])

    print("start")
    while time.time() < stop_time:
        for i in range(num_sensors):
            pi.write(pins_light[i], 1)
            time.sleep(1)
            timepoint, intensity_datapoint = determine_intensity_single_channel(pi,
                                                                                pins_light[i], i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels[i])
            intensity[i].append(intensity_datapoint)
            timepoints[i].append(timepoint - start_time)
            print(intensity)
    print("test completed, processing results")
    for i in range(num_sensors):  # check this
        intensity_single_channel = intensity[i]
        #absorbance[i] = calculate_absorbance(intensity_single_channel)

write_absorbance_csv(timepoints[0:num_sensors], absorbance)

print("finished processing result")
