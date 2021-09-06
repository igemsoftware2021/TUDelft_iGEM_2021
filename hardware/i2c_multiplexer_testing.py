import time
from apds9930 import APDS9930
from apds9930.values import APDS9930_ATIME
import pigpio

# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()

i2c_multiplexer_adress = 0x70
i2c_multiplexer_handle = pi.i2c_open(1, i2c_multiplexer_adress)

pi.i2c_write_device(i2c_multiplexer_handle, [0x80 | 0x01])


time.sleep(0.12)

i2c_sensor_handle = APDS9930(1)
i2c_sensor_handle.enable_ambient_light_sensor()
#i2c_sensor_handle.ambient_light_gain = 1
i2c_sensor_handle.write_byte_data(
    APDS9930_ATIME, 0xc0)  # Integration time 175 ms


print(type(i2c_sensor_handle))
stop = time.time() + 60
while time.time() < stop:
    i2c_sensor_handle.ch0_light
    i2c_sensor_handle.ch1_light
    result = i2c_sensor_handle.ambient_light
    result_0 = i2c_sensor_handle.ch0_light
    result_1 = i2c_sensor_handle.ch1_light
    print("Lux", result)
    print("CH0", result_0)
    print("CH1", result_1)
    time.sleep(0.3)
    print(time.time() - stop + 60)
    if time.time() > stop - 55:
        pi.i2c_write_device(i2c_multiplexer_handle, [0x80 | 0x02])

# TODO van te voren testen zit je te cappen op CH0 met een blanco, overweeg anders andere channel

i2c_sensor_handle.close()
pi.close()
