import time
import numpy as np
import pigpio
import csv


def i2c_multiplexer_select_channel(pi, i2c_multiplexer_handle, channel_number):
    """Function to select a channel on the TCA-9548 I2C multiplexer.
    """
    channel_number_base_2 = 2 ** channel_number  # bit on position determines channel on
    pi.i2c_write_device(i2c_multiplexer_handle,
                        [(channel_number_base_2)])
    print(channel_number)


def i2c_write_to_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, reg, data):
    """Function to write a data to a spefic registor of the ALS sensors that are connected to the multiplexer.
    """
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(pi,
                                       i2c_multiplexer_handle, channel_number)
        i2c_sensor_handle.write_byte_data(
            reg, data)


def i2c_set_gain_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, i2c_sensor_gain):
    """Function to change the gain of the ALS sensors that are connected to the multiplexer.
    The variable i2c_sensor_gain can take the values 0 (gain = 1), 1 (gain = 8), 2 (gain = 16), and 3 (gain = 130).
    """
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(pi,
                                       i2c_multiplexer_handle, channel_number)
        i2c_sensor_handle.ambient_light_gain = i2c_sensor_gain


def i2c_change_gain_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, gain):
    """Function to change the gain of the ALS sensors that are connected to the multiplexer.
    variable i2c_sensor_gain can take the values 0 (gain = 1), 1 (gain = 8), 2 (gain = 16), and 3 (gain = 130).
    """
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(
            i2c_multiplexer_select_channel, channel_number)
        i2c_sensor_handle.ambient_light_gain = gain


def i2c_activate_als_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers):
    """Function to initialize the ALS sensing of the APDS9930 sensors.
    """
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(pi,
                                       i2c_multiplexer_handle, channel_number)
        i2c_sensor_handle.enable_ambient_light_sensor()


def determine_intensity_single_channel(pi, pin_light, i2c_multiplexer_handle, i2c_sensor_handle, channel_number):
    """Function to measure the intensity of a single ambient light sensor.
    """
    pi.write(pin_light, 1)
    i2c_multiplexer_select_channel(pi,
                                   i2c_multiplexer_handle, channel_number)
    intensity = i2c_sensor_handle.ch0_light
    timepoint = time.time()
    time.sleep(0.25)
    pi.write(pin_light, 0)
    return timepoint, intensity


def calculate_absorbance(intensity):
    """Function to calculate the absorbance using the Beer-Lambert law.
    """
    blanco = intensity[0]
    intensity_np_arr = np.array(intensity)
    transmittance = intensity_np_arr / blanco
    absorbance = - np.log10(transmittance)
    return absorbance


def read_mcp3008(pi, adc, channel):
    """Function to calculate the absorbance using the Beer-Lambert law.
    Adapted from https://forums.raspberrypi.com/viewtopic.php?t=249872.
    """
    count, data = pi.spi_xfer(adc, [1, (8 + channel) << 4, 0])
    value = ((data[1] << 8) | data[2]) & 0x3FF
    return value


def read_mcp3008_median(pi, adc, channel, v_ref, gain, duration=0.17, interval=0.005):
    """Function to read-out the measurement of the temperature sensor by communicating with the ADC-converter (MCP3008) using the SPI-protocol.
    """
    samples = int(np.floor(duration / interval))
    adc_values = np.zeros(samples, dtype=np.float32)
    for i in range(samples):
        adc_values[i] = read_mcp3008(pi, adc, channel)
        time.sleep(interval)
    adc_values_median = np.median(adc_values)
    adc_voltage = adc_values_median * v_ref / 1024
    adc_temperature = adc_voltage * 100 / gain
    return adc_temperature


def pid_controller_calculator(pid_parameters, error, duration=0.17):
    """Function to calculate the output of a Proportional - Integral - Derivative - Controller (PID-controller).
    """
    k_p = pid_parameters[0]
    k_i = pid_parameters[1]
    k_d = pid_parameters[2]
    error = np.array(error)
    proportional_term = k_p * error[-1]
    integrational_term = k_i * duration * np.sum(error)
    differential_term = k_d * (error[-1] - error[-2]) / duration
    output = np.sqrt(
        (proportional_term + integrational_term + differential_term) * 22)
    return output


def check_heated_up(error):
    """Function to check if the test box is stably heated to the desired temperature.
    """
    error = np.array(error)
    error = np.abs(error)
    if len(error) < 30 or sum(error[-30:-1]) > 15:
        return False
    else:
        return True


def pre_heater(pi, adc, channel, pin_heating, pid_parameters, v_ref, gain, duration=0.17, interval=0.005):
    """Function to pre-heat the test box to the desired temperature before the test is begun.
    """
    temperature_desired = pid_parameters[3]
    pwm_freq = pid_parameters[4]
    temperature_error = []
    temperature_time = []
    heated_up = False

    T_c = 0.001  # Lowest high or low time [ms] to protect the hardware
    duty_cycle_lower_bound = T_c * pwm_freq * 10 ** 6
    duty_cycle_upper_bound = (1 - T_c * pwm_freq) * 10**6

    while heated_up == False:
        # Measure temperature error
        temperature_median = read_mcp3008_median(pi, adc, channel, v_ref, gain,
                                                 duration=duration, interval=interval)
        temperature_time.append(time.time())
        temperature_median_error = temperature_desired - temperature_median
        temperature_error.append(temperature_median_error)

        # PID controller
        duty_cycle = pid_controller_calculator(pid_parameters, temperature_error,
                                               duration=duration)
        duty_cycle = int(duty_cycle)  # PWM will have freq
        if duty_cycle < duty_cycle_lower_bound:
            duty_cycle = 0
        if duty_cycle > duty_cycle_upper_bound:
            duty_cycle = 1*10**6
        print(duty_cycle)
        pi.hardware_PWM(pin_heating, pwm_freq, duty_cycle)

        # Register that and display message if the box is stable at desired temperature
        heated_up = check_heated_up(temperature_error)
        if heated_up == True:
            # Turn of the PWM output signal
            pi.hardware_PWM(pin_heating, pwm_freq, 0)
            print("pre-heating has been completed")
            return temperature_time, temperature_error


def initialize_light_pins(pi, pins):
    """Function to set the output mode and the pull up resistor of the GPIO pins used to control the LEDs.
    """
    for pin in pins:
        pi.set_mode(pin, pigpio.OUTPUT)
        pi.set_pull_up_down(pin, pigpio.PUD_DOWN)


def initialize_heating_pin(pi, pin):
    """Function to initialize the pull up resistor on the GPIO to switch the heating elements on and off.
    """
    pi.set_pull_up_down(pin, pigpio.PUD_DOWN)


def temperature_controller(pi, adc, channel, pin_heating, duty_cycle_lower_bound, duty_cycle_upper_bound, pid_parameters, v_ref, gain, temperature_error, duration=0.17, interval=0.005):
    """Function to keep the temperature at the desired temperature during the test.
    """
    temperature_desired = pid_parameters[3]
    pwm_freq = pid_parameters[4]

    temperature_median = read_mcp3008_median(pi, adc, channel, v_ref, gain,
                                             duration=duration, interval=interval)
    temperature_median_error = temperature_desired - temperature_median
    temperature_error.append(temperature_median_error)

    # PID controller
    duty_cycle = pid_controller_calculator(pid_parameters, temperature_error,
                                           duration=duration)
    duty_cycle = int(duty_cycle)  # PWM will have freq
    if duty_cycle < duty_cycle_lower_bound:
        duty_cycle = 0
    if duty_cycle > duty_cycle_upper_bound:
        duty_cycle = duty_cycle_upper_bound
        print(duty_cycle)
    pi.hardware_PWM(pin_heating, pwm_freq, duty_cycle)
    return temperature_error
