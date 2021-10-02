"""
This module create pseudo random binary sequence (prbs) with the polynomial  x^7 + x^6 +1.
"""
import numpy as np


def bin_array(num, m):
    """Convert a positive integer num into an m-bit bit vector"""
    return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8)


start = np.uint8(0X02)
a = start
print(start)

prbs_7 = []
i = 0
while 1:
    newbit = (((a >> 6) ^ (a >> 5)) & 1)
    a = ((a << 1) | newbit) & 0x7f
    print(a)
    i += 1
    if a == start:
        print(i)
        break

print(i)
