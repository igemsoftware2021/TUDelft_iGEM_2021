import time
import pigpio

# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()

# Initialize the light control pins
#pin_light_1 = 17

pin_transistor_1 = 27

#pi.set_mode(pin_light_1, pigpio.OUTPUT)
pi.set_mode(pin_transistor_1, pigpio.OUTPUT)


for i in range(3):
    # Set the light pin as high
    pi.write(pin_transistor_1, 1)

    # Some delay to see the effect
    time.sleep(10)

    # Turn off the light again
    pi.write(pin_transistor_1, 0)

    # Some delay to see the effect
    time.sleep(1)

# Close down the pi
pi.stop()
