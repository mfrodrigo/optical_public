"""

"""
from numpy import repeat
from numpy.random import randint


def build_square_wave(num_samplesperbit, size=None, input_seq=None):
    if input_seq is None:
       input_seq = randint(2, size=size).reshape(size, 1)
    return repeat(input_seq, num_samplesperbit, axis=0)
