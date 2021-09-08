from hardware.helpers import determine_intensity_all_channels
import time
import multiprocessing
import numpy as np
import pigpio
import helpers

###############################################################################
# Set-up

# Constants for test
total_time = 10

# Control pins
pin_heating = 27
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
pid_parameters = [k_p, k_i, k_d, temperature_desired]
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
helpers.initialize_light_pins(pi, pins_light)
helpers.initialize_heating_pin(pi, pin_heating, pwm_freq)


# Initialize SPI connection
adc_handle = pi.spi_open(spi_channel, spi_baud, spi_flags)

# Initialize I2C connection
i2c_multiplexer_handle = pi.i2c_open(i2c_bus, i2c_multiplexer_adress)
helpers.i2c_multiplexer_select_channel(pi,
                                       i2c_multiplexer_handle, i2c_sensor_1_channel)  # Switch to a sensor channel
i2c_sensor_handle = APDS9930(i2c_bus)  # Create 1 sensor handle

# Activate ALS, and set gain and integration time
helpers.i2c_activate_als_all_sensors(
    pi, i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels)
helpers.i2c_write_to_all_sensors(i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels,
                                 APDS9930_ATIME, i2c_sensor_int_time)
# helpers.i2c_change_gain_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, gain)

# TODO test gain

###############################################################################


print("Would you like to pre-heat the box? (True or False)")
pre_heating = input()
if pre_heating == True:
    helpers.pre_heater(pi, adc_handle, spi_channel, pin_heating,
                       pid_parameters, v_ref, gain)
    print("Pre-heating has finished.")


print("Would you like to run a test? (True or False)")
test = input()
if test == True:
    p_heating = multiprocessing.Process(target=helpers.temperature_controller, args=(pi,
                                                                                     adc_handle, spi_channel, total_time, pin_heating, pid_parameters, v_ref, gain))
    p_light = multiprocessing.Process(target=helpers.determine_intensity_over_time, args=(pi,
                                                                                          pins_light, i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels, total_time))
    p_heating.start()
    p_light.start()
    p_heating.join()
    p_light.join()
print("The test has been completed.")
