from copy import error
import time
import multiprocessing
import numpy as np
import pigpio


def i2c_multiplexer_select_channel(pi, i2c_multiplexer_handle, channel_number):
    channel_decimal = 2 ** channel_number
    channel_hex = hex(channel_decimal)
    print(type(channel_hex))
    print(type(0x80))
    pi.i2c_write_device(i2c_multiplexer_handle, [0x80 | channel_hex])
    time.sleep(0.1)


def i2c_write_to_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, reg, data):
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(pi,
                                       i2c_multiplexer_handle, channel_number)
        i2c_sensor_handle.write_byte_data(
            reg, data)


def i2c_change_gain_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, gain):
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(
            i2c_multiplexer_select_channel, channel_number)
        i2c_sensor_handle.ambient_light_gain = gain


def i2c_activate_als_all_sensors(pi, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers):
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(
            i2c_multiplexer_select_channel, channel_number)
        i2c_sensor_handle.enable_ambient_light_sensor()


def determine_intensity_single_channel(pi, pin_light, i2c_multiplexer_handle, i2c_sensor_handle, channel_number):
    pi.write(pin_light, 1)
    i2c_multiplexer_select_channel(pi,
                                   i2c_multiplexer_handle, channel_number)
    intensity = i2c_sensor_handle.ch0_light
    timepoint = time.time()
    return timepoint, intensity


def determine_intensity_over_time(pi, pins_light, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, total_time):
    num_sensors = len(channel_numbers)
    start_time = time.time()
    stop_time = start_time() + total_time
    intensity = [[]] * num_sensors
    timepoints = []
    for i in range(num_sensors):
        timepoints.append([])
    while time.time() < stop_time:
        for i in range(len(num_sensors)):
            pi.write(pins_light[i])
            i2c_multiplexer_select_channel(pi,
                                           i2c_multiplexer_handle, channel_numbers[i])
            timepoint, intensity_datapoint = determine_intensity_single_channel(pi,
                                                                                pins_light[i], i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers[i])
            intensity[i].append(intensity_datapoint)
            timepoints[i].append(timepoint)
    for i in range(len(num_sensors)):
        intensity_single_channel = intensity[i]
        timepoints_single_channel = timepoints[i]
        absorbance_single_channel = calculate_absorbance(
            intensity_single_channel)
        # TODO write to file with something to denote well
    return timepoints, intensity


def calculate_absorbance(intensity):
    blanc = intensity[0]
    intensity_np_arr = np.array(intensity)
    ratio = intensity_np_arr / blanc
    absorbance = np.log10(ratio)
    return absorbance


def read_mcp3008(pi, adc, channel):  # TODO link to stackoverflow
    count, data = pi.spi_xfer(adc, [1, (8 + channel) << 4, 0])
    value = ((data[1] << 8) | data[2]) & 0x3FF
    return value


def read_mcp3008_median(pi, adc, channel, v_ref, gain, duration=0.25, interval=0.025):
    samples = duration / interval
    adc_values = np.zeros(samples, dtype=np.float32)
    for i in range(samples):
        adc_values[i] = read_mcp3008(pi, adc, channel)
        time.sleep(interval)
    adc_values_median = np.median(adc_values)
    adc_voltage = adc_values_median * v_ref / 1024
    adc_temperature = adc_voltage * 100 / gain
    return adc_temperature


def pid_controller_calculator(pid_parameters, error, duration):
    k_p = pid_parameters[0]
    k_i = pid_parameters[1]
    k_d = pid_parameters[2]
    error = np.array(error)
    proportional_term = k_p * error[-1]
    integrational_term = k_i * duration * np.sum(error)
    differential_term = k_d * 0  # TODO adapt
    output = proportional_term + integrational_term + differential_term
    return output


def check_heated_up(error):
    error = np.array(error)
    error = np.abs(error)
    if len(error) < 30 or sum(error[-30:-1]) > 15:
        return False
    else:
        return True


def pre_heater(pi, adc, channel, pin_heating, pid_parameters, v_ref, gain, duration=0.25, interval=0.005):
    temperature_desired = pid_parameters[4]
    temperature_error = []
    heated_up = False

    while heated_up == False:
        # Measure temperature error
        temperature_median = read_mcp3008_median(pi, adc, channel, v_ref, gain,
                                                 duration=duration, interval=interval)
        temperature_median_error = temperature_desired - temperature_median
        temperature_error.append(temperature_median_error)

        # PID controller
        duty_cycle = pid_controller_calculator(pid_parameters, temperature_error,
                                               duration=0.25)
        duty_cycle = int(duty_cycle)
        pi.set_PWM_dutycycle(pin_heating, duty_cycle)

        # Register that and display message if the box is stable at desired temperature
        heated_up = check_heated_up(temperature_error)
        if heated_up == True:
            return


def temperature_controller(pi, adc, channel, total_time, pin_heating, pid_parameters, v_ref, gain, duration=0.25, interval=0.005):
    stop_time = time.time() + total_time
    temperature_desired = pid_parameters[4]
    temperature_error = []

    while time.time() < stop_time:
        # Measure temperature error
        temperature_median = read_mcp3008_median(pi, adc, channel, v_ref, gain,
                                                 duration=duration, interval=interval)
        temperature_median_error = temperature_desired - temperature_median
        temperature_error.append(temperature_median_error)

        # PID controller
        duty_cycle = pid_controller_calculator(pid_parameters, temperature_error,
                                               duration=0.25)
        duty_cycle = int(duty_cycle)
        pi.set_PWM_dutycycle(pin_heating, duty_cycle)
    pi.set_PWM_dutycycle(pin_heating, 0)  # Turn of the PWM output signal
    return error


def initialize_light_pins(pi, pins):
    for pin in pins:
        pi.set_mode(pin, pigpio.OUTPUT)
        pi.set_pull_up_down(pin, pigpio.PUD_DOWN)


def initialize_heating_pin(pi, pin, pwm_freq):
    pi.set_pull_up_down(pin, pigpio.PUD_DOWN)
    pi.set_PWM_frequency(pin, pwm_freq)
