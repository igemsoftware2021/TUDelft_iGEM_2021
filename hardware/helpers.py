import time
import multiprocessing
import numpy as np
import pigpio


def i2c_multiplexer_select_channel(i2c_multiplexer_handle, channel_number):
    channel_decimal = 2 ** channel_number
    channel_hex = hex(channel_decimal)
    pi.i2c_write_device(i2c_multiplexer_handle, [0x80 | channel_hex])
    time.sleep(0.1)


def i2c_write_to_all_sensors(i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, reg, data):
    for channel_number in channel_numbers:
        i2c_multiplexer_select_channel(
            i2c_multiplexer_handle, channel_number)
        i2c_sensor_handle.write_byte_data(
            reg, data)


def determine_absorbance(pin_light, i2c_multiplexer_handle, i2c_sensor_handle, channel_number):
    pi.write(pin_light, 1)
    i2c_multiplexer_select_channel(
        i2c_multiplexer_handle, channel_number)
    intensity = i2c_sensor_handle.ch0_light
    timepoint = time.time()
    return timepoint, intensity


def determine_absorbance_over_time(pins_light, i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers, total_time):
    num_sensors = len(channel_numbers)
    start_time = time.time()
    stop_time = start_time() + total_time
    while time.time() < stop_time:
        for i in range(len(num_sensors)):
            pi.write(pins_light[i])
            i2c_multiplexer_select_channel(
                i2c_multiplexer_handle, channel_numbers[i])
            timepoint, intensity = determine_absorbance(
                pins_light[i], i2c_multiplexer_handle, i2c_sensor_handle, channel_numbers[i])
            # TODO store these in a list?
