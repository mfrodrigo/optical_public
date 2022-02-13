"""
Filter Design.
"""
from numpy import ndarray, pi
from scipy.signal import butter, filtfilt


def butter_filter(N, cutoff,  fs, sig) -> ndarray:
    """

    """
    wn = cutoff/(fs/2)
    b, a = butter(N, Wn=wn, btype='lowpass')
    return filtfilt(b, a, sig)
