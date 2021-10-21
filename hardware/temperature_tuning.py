import time
import pigpio
import numpy as np
import helpers
import matplotlib.pyplot as plt
from file_helpers import write_temperature_csv

# Heating pin
pin_heating = 13

# Constants for heating
k_u = 1000000  # Ultimate gain, used for Ziegler-Nichols tuning of the PID controller
temperature_desired = 37  # Desired temperature in degrees Celsius
pwm_freq = 10
controller_parameters = [k_u, temperature_desired, pwm_freq]
v_ref = 3.3
gain = 2.96

# SPI
spi_channel = 0  # SPI channel of the Raspberry Pi that is connected to the ADC
spi_baud = 50000  # Speed of the serial communication in bits/s
spi_flags = 0


# Initialize Raspberry Pi
pi = pigpio.pi()

# Initialize SPI connection
adc_handle = pi.spi_open(spi_channel, spi_baud, spi_flags)
# Initialize heating pin
helpers.initialize_heating_pin(pi, pin_heating)

pi.set_mode(pin_heating, pigpio.OUTPUT)
pi.write(pin_heating, 1)
print("licht")
time.sleep(10)

temperature_error = []
temperature_time = []

T_c = 0.001  # Lowest high or low time [ms] to protect the hardware
duty_cycle_lower_bound = T_c * pwm_freq * 10 ** 6
duty_cycle_upper_bound = (1 - T_c * pwm_freq) * 10**6
duration = 0.10
interval = 0.005
pi.hardware_PWM(pin_heating, pwm_freq, 10**6)
print("nu")
time.sleep(10)
start = time.time()
print("hoi")
print(start)
print(start+120)
while time.time() < start + 120:
    print("test")
    # Measure temperature error
    temperature_median = helpers.read_mcp3008_median(pi, adc_handle, spi_channel, v_ref, gain,
                                                     duration=duration, interval=interval)
    temperature_time.append(time.time())
    temperature_median_error_timepoint = temperature_desired - temperature_median
    temperature_error.append(temperature_median_error_timepoint)
    print(temperature_median)

    duty_cycle = np.sqrt(k_u * temperature_error[-1]/22)
    duty_cycle = int(duty_cycle)
    if duty_cycle < duty_cycle_lower_bound:
        duty_cycle = 0
    if duty_cycle > duty_cycle_upper_bound:
        duty_cycle = 1*10**6
    pi.hardware_PWM(pin_heating, pwm_freq, duty_cycle)


path = "hardware/data_tuning" + f"temperature-{time.time()}.csv"
write_temperature_csv(temperature_time, temperature_error, path)

if not pi.connected:
    exit()


# Close down the pi
pi.stop()
