"""
Filter Design.
"""
import numpy as np
from numpy import ndarray, pi, abs, log10
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft

from output.plotter import Plotter


def transfer_function(frequencies: ndarray) -> ndarray:
    """
    Calculate the H(f) for a 4 order
    butterworth filter and cutoff frequency equal to 1 GHz.
    Args:
        frequencies: (ndarray) frequency vector

    Returns:
        H: (ndarray) Response in a certain region of filter frequencies.
    """
    s = 1j * 2 * pi * frequencies
    H = 1.559e39 / (s ** 4 + 1.642e10 * s ** 3 + 1.348e20 * s ** 2 + 6.482e29 * s + 1.559e39)
    return H
