import time
import numpy as np
import pigpio
import csv


def i2c_multiplexer_select_channel(pi, i2c_multiplexer_handle, channel_number):
    channel_base_2_number = 2 ** channel_number
    pi.i2c_write_device(i2c_multiplexer_handle,
                        [0x80 | channel_base_2_number])


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
        i2c_multiplexer_select_channel(pi,
                                       i2c_multiplexer_handle, channel_number)
        i2c_sensor_handle.enable_ambient_light_sensor()


def determine_intensity_single_channel(pi, pin_light, i2c_multiplexer_handle, i2c_sensor_handle, channel_number):
    pi.write(pin_light, 1)
    i2c_multiplexer_select_channel(pi,
                                   i2c_multiplexer_handle, channel_number)
    intensity = i2c_sensor_handle.ch0_light
    timepoint = time.time()
    pi.write(pin_light, 0)
    return timepoint, intensity


def calculate_absorbance(intensity):
    blanco = intensity[0]
    intensity_np_arr = np.array(intensity)
    transmittance = intensity_np_arr / blanco
    absorbance = - np.log10(transmittance)
    return absorbance


def read_mcp3008(pi, adc, channel):  # TODO link to stackoverflow
    count, data = pi.spi_xfer(adc, [1, (8 + channel) << 4, 0])
    value = ((data[1] << 8) | data[2]) & 0x3FF
    return value


def read_mcp3008_median(pi, adc, channel, v_ref, gain, duration=0.17, interval=0.005):
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


def pre_heater(pi, adc, channel, pin_heating, pid_parameters, v_ref, gain, duration=0.17, interval=0.005):
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
            duty_cycle = duty_cycle_upper_bound
        print(duty_cycle)
        pi.hardware_PWM(pin_heating, pwm_freq, duty_cycle)

        # Register that and display message if the box is stable at desired temperature
        heated_up = check_heated_up(temperature_error)
        if heated_up == True:
            # Turn of the PWM output signal
            pi.hardware_PWM(pin_heating, pwm_freq, 0)
            return temperature_time, temperature_error


def initialize_light_pins(pi, pins):
    for pin in pins:
        pi.set_mode(pin, pigpio.OUTPUT)
        pi.set_pull_up_down(pin, pigpio.PUD_DOWN)


def initialize_heating_pin(pi, pin):
    pi.set_pull_up_down(pin, pigpio.PUD_DOWN)


def main_measurement(pi, pins_light, pin_heating, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, adc_handle, adc_channel, total_time, pid_parameters, v_ref, gain, duration=0.17, interval=0.005):
    num_sensors = len(channel_numbers)
    pwm_freq = pid_parameters[4]
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
            temperature_controller(pi, adc_handle, adc_channel, pin_heating, duty_cycle_lower_bound, duty_cycle_upper_bound,
                                   pid_parameters, v_ref, gain, temperature_error, duration=duration, interval=interval)
            temperature_timepoint = time.time()
            timepoint, intensity_datapoint = determine_intensity_single_channel(pi,
                                                                                pins_light[i], i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers[i])
            intensity[i].append(intensity_datapoint)
            timepoints[i].append(timepoint - start_time)
            timepoints[-1].append(temperature_timepoint - start_time)

    pi.hardware_PWM(pin_heating, pwm_freq, 0)  # Turn of the PWM output signal

    for i in range(num_sensors):
        intensity_single_channel = intensity[i]
        absorbance[i] = calculate_absorbance(intensity_single_channel)
        print(intensity_single_channel)
    return timepoints, absorbance, temperature_error


def temperature_controller(pi, adc, channel, pin_heating, duty_cycle_lower_bound, duty_cycle_upper_bound, pid_parameters, v_ref, gain, temperature_error, duration=0.17, interval=0.005):
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


def save_as_csv(timepoints, absorbance, path, name):
    timepoints = timepoints
    absorbance = absorbance
    num_rows = len(timepoints)
    print(num_rows)
    with open(path + name + ".csv", 'w') as new_file:
        csv_writer = csv.writer(new_file)
        for i in range(num_rows):
            row = [str(timepoints[i]), str(absorbance[i])]
            print(row)
            csv_writer.writerow(row)


def read_from_csv(path, name):
    timepoints = []
    absorbance = []
    with open(path + name + ".csv", "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for line in csv_reader:
            print(line)
            timepoints.append(np.float32(line[0]))
            absorbance.append(np.float32(line[1]))
    timepoints = np.array(timepoints)
    absorbance = np.array(absorbance)
    return timepoints, absorbance
