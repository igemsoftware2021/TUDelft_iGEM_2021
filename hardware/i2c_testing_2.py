import time
import pigpio
from smbus2 import SMBus

# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()

bus_1 = SMBus(1)

i2c_addr = 0x39

# Constants
ATIME = 0xf6
WTIME = 0xff
PTIME = 0xff
PPULSE = 0x01


bus_1.write_byte_data(i2c_addr, 0x00, 0x00)
#bus_1.write_byte_data(i2c_addr,  0x01, ATIME)
bus_1.write_byte_data(i2c_addr,  0x02, PTIME)
bus_1.write_byte_data(i2c_addr,  0x03, WTIME)
bus_1.write_byte_data(i2c_addr,  0x0E, PPULSE)

PDRIVE = 0x00
PDIODE = 0x20
PGAIN = 0x00
AGAIN = 0x00

bus_1.write_byte_data(i2c_addr,  0x0F, PDRIVE | PDIODE | PGAIN | AGAIN)

WEN = 0x08
PEN = 0x04
AEN = 0x02
PON = 0x01

bus_1.write_byte_data(i2c_addr,  0x00, WEN | PEN | AEN | PON)

time.sleep(0.12)


print(bus_1.read_i2c_block_data(i2c_addr, 0x01, 10))

bus_1.close()


# Close down the pi
pi.stop()
