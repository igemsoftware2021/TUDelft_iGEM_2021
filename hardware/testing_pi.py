import time
import pigpio

# Initialize Raspberry Pi
pi = pigpio.pi()

# Initialize the light control pins
pin_light_1 = 17

pi.set_mode(pin_light_1, pigpio.OUTPUT)


for i in range(10):
    # Set the light pin as high
    pi.write(pin_light_1, 1)

    # Some delay to see the effect
    time.sleep(0.1)

    # Turn off the light again
    pi.write(pin_light_1, 0)

    # Some delay to see the effect
    time.sleep(0.1)

# Close down the pi
pi.stop()
