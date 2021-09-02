import time
import pigpio
import numpy as np

# Initialize Raspberry Pi
pi = pigpio.pi()

# Initialize SPI with ADC
# spi_flags = [0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
#   flags     [b, b, b, b, b, b, R, T, n, n, n, n, W, A, u2,u1,u0,p2,p1,p0, m, m]
spi_flags = 688576
adc = pi.spi_open(0, 50000, spi_flags)
print("Opening succesfull")
pi.spi_write(adc, [128])
print("Writing succesfull")

n = 100
v_ref = 5
n_measured, adc_bytes = pi.spi_read(adc, n)

pi.spi_close(adc)

# Close down the pi
pi.stop()

print(adc_bytes)
adc_integers = np.zeros(n)
if n_measured == 100:
    print("Reading succesfull")
    for i in range(n):
        adc_integers[i] = int.from_bytes(adc_bytes[i], "big")
        adc_voltage = adc_integers * v_ref / 1024
        adc_temperature = adc_voltage * 100
        adc_temperature_mean = np.mean(adc_temperature)
        print(adc_temperature_mean)
else:
    print("Error: Data acquisition did not succeed.")
