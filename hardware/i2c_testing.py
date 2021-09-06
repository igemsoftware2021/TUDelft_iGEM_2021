import time
import pigpio


# Initialize Raspberry Pi
pi = pigpio.pi()

if not pi.connected:
    exit()


def adps9930_read_als(pi, handle):
    num_low = -1
    num_high = -1
    while num_low < 0:
        num_low, low = pi.i2c_read_i2c_block_data(handle, 0x14, 2)
    while num_high < 0:
        num_high, high = pi.i2c_read_i2c_block_data(handle, 0x15, 2)
    print(low)
    print(low[1])
    print(high)
    print(high[1])
    bytevalue = high[1] << 8  # First shifted 8 bits
    print(bytevalue)
    bytevalue = bytevalue | low[1]  # Then bitwise with low
    print(bytevalue)
    # Because we are working in hexidecimals, base=16
    value = int(bytevalue, base=16)
    return value


i2c_handle = pi.i2c_open(1, 0x39)

# print(f"0x01: {pi.i2c_read_i2c_block_data(i2c_handle, 0x01, 4)}")
# print(f"0xA0 | 0x01: {pi.i2c_read_i2c_block_data(i2c_handle, 0xA0 | 0x01, 4)}")
# # print(pi.i2c_read_i2c_block_data(i2c_handle, 0x01, 10))


# pi.i2c_write_byte_data(i2c_handle, 0x00, 0x00)
# #pi.i2c_write_device(i2c_handle, [0x80 | 0x01, ATIME])

# print(f"0x01: {pi.i2c_read_i2c_block_data(i2c_handle, 0x01, 4)}")
# print(f"0xA0 | 0x01: {pi.i2c_read_i2c_block_data(i2c_handle, 0xA0 | 0x01, 4)}")
# print("Now writing")
# pi.i2c_write_i2c_block_data(i2c_handle, 0x80 | 0x01, [ATIME])
# print(f"0x01: {pi.i2c_read_i2c_block_data(i2c_handle, 0x01, 4)}")
# print(f"0xA0 | 0x01: {pi.i2c_read_i2c_block_data(i2c_handle, 0xA0 | 0x01, 4)}")

# print("second thingie")
# result1 = pi.i2c_read_device(i2c_handle, 10)
# print(result1)
# pi.i2c_write_device(i2c_handle, [0x80 | 0x01, ATIME])
# print("writing")
# result2 = pi.i2c_read_device(i2c_handle, 10)
# print(result2)

# print("second 2.0 thingie")
# result1 = pi.i2c_read_device(i2c_handle, 1)
# print(result1)
# pi.i2c_write_device(i2c_handle, [0x80 | 0x02, WTIME])
# print("writing")
# result2 = pi.i2c_read_device(i2c_handle, 1)
# print(result2)

# print("third thingie")
# result1 = pi.i2c_read_device(i2c_handle, 10)
# print(result1)
# #pi.i2c_write_device(i2c_handle, [0x14])
# print("writing")
# result2 = pi.i2c_read_byte_data(i2c_handle, 0xA0 | 0x14)
# result3 = pi.i2c_read_byte_data(i2c_handle, 0xA0 | 0x15)
# print(result2)
# print(result3)

# pi.i2c_write_byte_data(i2c_handle, 0x01, ATIME)
# pi.i2c_write_byte_data(i2c_handle, 0x02, PTIME)
# pi.i2c_write_byte_data(i2c_handle, 0x03, WTIME)
# pi.i2c_write_byte_data(i2c_handle, 0x80 | 0x0E, PPULSE)
# result4 = pi.i2c_read_device(i2c_handle, 1)
# print("result4: ")
# print(result4)

# Constants
ATIME = 0xf6
WTIME = 0xb6
PTIME = 0xff
PPULSE = 0x01

pi.i2c_write_device(i2c_handle, [0x80 | 0x01, ATIME])
pi.i2c_write_device(i2c_handle, [0x80 | 0x02, WTIME])
pi.i2c_write_device(i2c_handle, [0x80 | 0x03, PTIME])
pi.i2c_write_device(i2c_handle, [0x80 | 0x0E, PPULSE])


PDRIVE = 0x00
PDIODE = 0x20
PGAIN = 0x00
AGAIN = 0x00

pi.i2c_write_device(i2c_handle, [0x80 | 0x0F, PDRIVE | PDIODE | PGAIN | AGAIN])

WEN = 0x08
PEN = 0x04
AEN = 0x02
PON = 0x01

pi.i2c_write_device(i2c_handle, [0x80 | 0x00, WEN | PEN | AEN | PON])

time.sleep(0.012)
stop = time.time() + 20

while time.time() < stop:
    pi.i2c_write_device(i2c_handle, [0xA0 | 0x14])
    num, result_low = pi.i2c_read_device(i2c_handle, 5)

    pi.i2c_write_device(i2c_handle, [0xA0 | 0x15])
    num, result_high = pi.i2c_read_device(i2c_handle, 5)

    print(result_low, result_high)
    time.sleep(0.3)

#print(pi.i2c_read_i2c_block_data(i2c_handle, 0x14, 10))
#print(pi.i2c_read_byte_data(i2c_handle, 0x14))
stop = time.time() + 15
exit()

while time.time() < stop:
    time.sleep(0.1)
    als_value = adps9930_read_als(pi, i2c_handle)
    print(als_value)


pi.i2c_close(i2c_handle)


# Close down the pi
pi.stop()
