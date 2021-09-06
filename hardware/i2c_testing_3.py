import time
from apds9930 import APDS9930

i2c_handle = APDS9930(1)
i2c_handle.enable_ambient_light_sensor()

print(type(i2c_handle))
stop = time.time() + 60
while time.time() < stop:
    i2c_handle.ch0_light
    i2c_handle.ch1_light
    result = i2c_handle.ambient_light
    print(result)
    time.sleep(0.1)

i2c_handle.close()
