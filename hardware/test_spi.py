import time
import pigpio
import numpy as np


def read_mcp3008(adc, channel):
    count, data = pi.spi_xfer(adc, [1, (8 + channel) << 4, 0])
    value = ((data[1] << 8) | data[2]) & 0x3FF
    return value


# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()


v_ref = 3.3  # Reference voltage of the ADC-converter.
gain = 2.96  # Gain of the amplifier connected to the temperature sensor

adc = pi.spi_open(0, 50000, 0)  # CE0 (indicates to which slave it's talking)

stop = time.time() + 5
while time.time() < stop:
    adc_value = read_mcp3008(adc, 0)
    print(adc_value)
    adc_voltage = adc_value * v_ref / 1024
    adc_temperature = adc_voltage * 100 / gain
    print(adc_temperature)


pi.spi_close(adc)

# Close down the pi
pi.stop()
