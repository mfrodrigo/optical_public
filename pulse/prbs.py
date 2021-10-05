"""
This module create pseudo random binary sequence (prbs) with the polynomial  x^7 + x^6 +1.
"""
import numpy as np


def bin_array(num, m):
    """Convert a positive integer num into an m-bit bit vector"""
    return list(np.binary_repr(num).zfill(m))


def prbs_7(amount_of_values=128):
    start = np.uint8(0X02)
    a = start
    print(start)

    prbs_7 = []
    i = 0
    while 1:
        newbit = (((a >> 6) ^ (a >> 5)) & 1)
        a = ((a << 1) | newbit) & 0x7f
        prbs_7 = prbs_7 + bin_array(a, 7)
        i += 1
        if a == start:
            prbs_7 = np.array(prbs_7).astype(np.int8)
            break

    return prbs_7[:amount_of_values*7].reshape((len(prbs_7), 1))

