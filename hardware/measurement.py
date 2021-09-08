from hardware.helpers import i2c_multiplexer_select_channel
import time
import multiprocessing
import numpy as np
import pigpio
import helpers

###############################################################################
# Set-up

# Control pins
pin_heating = "TBD"
pin_light_1 = "TBD"
pin_light_2 = "TBD"
pin_light_3 = "TBD"
pin_light_4 = "TBD"
pins_light = [pin_light_1, pin_light_2, pin_light_3, pin_light_4]

# Constants for heating
k_u = "TBD"  # Ultimate gain, used for Ziegler-Nichols tuning of the PID controller
t_u = "TBD"  # Oscillation frequency associated with the ultimate gain
k_p = 0.20 * k_u        # Proportional gain PID controller
k_i = 0.50 * t_u        # Integral gain PID controller
k_d = 0.40 * k_u / t_u  # Derivative gain PID controller
pwm_freq = 10  # simulations suggest this is fine, see LT spice. Need to do this properly with a good toggle speed

# SPI
spi_channel = 0  # SPI channel of the Raspberry Pi that is connected to the ADC
spi_baud = 50000  # Speed of the serial communication in bits/s
spi_flags = 0

# I2C
i2c_bus = 1  # I2C bus of the Raspberry Pi that is connected to the multiplexer
i2c_multiplexer_adress = 0x70  # TCA9548
i2c_sensor_adress = 0x39  # APDS9930
i2c_sensor_int_time = 0xc0  # Integration time APDS9930. 0x0c <-> 175 ms
i2c_sensor_gain = "TBD"  # Gain of the APDS9930. #TODO implement gain can be changed
i2c_sensor_1_channel = 0  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_2_channel = 1  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_3_channel = 2  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_4_channel = 3  # Channel of the multiplexer to which the sensor is connected
i2c_sensor_channels = [i2c_sensor_1_channel, i2c_sensor_2_channel,
                       i2c_sensor_3_channel, i2c_sensor_4_channel]
###############################################################################


# Initialization

# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()
    # TODO display error

# Initialize pins
pi.set_mode(pin_light_1, pigpio.OUTPUT)
pi.set_mode(pin_light_2, pigpio.OUTPUT)
pi.set_mode(pin_light_3, pigpio.OUTPUT)
pi.set_mode(pin_light_4, pigpio.OUTPUT)

pi.set_pull_up_down(pin_light_1, pigpio.PUD_DOWN)
pi.set_pull_up_down(pin_light_2, pigpio.PUD_DOWN)
pi.set_pull_up_down(pin_light_3, pigpio.PUD_DOWN)
pi.set_pull_up_down(pin_light_4, pigpio.PUD_DOWN)

pi.set_pull_up_down(pin_heating, pigpio.PUD_DOWN)
pi.set_PWM_frequency(pin_heating, pwm_freq)


# Initialize SPI connection
adc_handle = pi.spi_open(spi_channel, spi_baud, spi_flags)

# Initialize I2C connection
i2c_multiplexer_handle = pi.i2c_open(i2c_bus, i2c_multiplexer_adress)
helpers.i2c_multiplexer_select_channel(
    i2c_multiplexer_handle, i2c_sensor_1_channel)  # Switch to a sensor channel
i2c_sensor_handle = APDS9930(i2c_bus)  # Create 1 sensor handle

for i in range(len(i2c_sensor_channels)):  # Activate ambient light sensing
    helpers.i2c_multiplexer_select_channel(
        i2c_multiplexer_select_channel, i2c_sensor_channels[i])
    i2c_sensor_handle.enable_ambient_light_sensor()

helpers.i2c_write_to_all_sensors(i2c_multiplexer_handle, i2c_sensor_handle, i2c_sensor_channels,
                                 APDS9930_ATIME, i2c_sensor_int_time)  # Set integration time


###############################################################################
