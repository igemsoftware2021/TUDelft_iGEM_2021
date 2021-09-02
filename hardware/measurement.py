import time
import multiprocessing
import numpy as np
import pigpio

# Initialize Raspberry Pi
pi = pigpio.pi()

# Initialize LED control pins
pin_light_1 = 17
pin_light_2 = 27
pin_light_3 = 22
pin_light_4 = 23

pi.set_mode(pin_light_1, pigpio.OUTPUT)
pi.set_mode(pin_light_2, pigpio.OUTPUT)
pi.set_mode(pin_light_3, pigpio.OUTPUT)
pi.set_mode(pin_light_4, pigpio.OUTPUT)

pi.set_pull_up_down(pin_light_1, pigpio.PUD_DOWN)
pi.set_pull_up_down(pin_light_2, pigpio.PUD_DOWN)
pi.set_pull_up_down(pin_light_3, pigpio.PUD_DOWN)
pi.set_pull_up_down(pin_light_4, pigpio.PUD_DOWN)
