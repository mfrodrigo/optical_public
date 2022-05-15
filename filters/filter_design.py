"""
Filter Design.
"""
from numpy import ndarray
from scipy.signal import butter, filtfilt, bessel


def butter_filter(N: int, cutoff: float, fs: float, sig: ndarray) -> ndarray:
    wn = cutoff / (fs / 2)
    b, a = butter(N, Wn=wn, btype='lowpass')
    return filtfilt(b, a, sig)


def bessel_filter(N: int, cutoff: float, fs: float, sig: ndarray) -> ndarray:
    wn = cutoff / (fs / 2)
    b, a = bessel(N, Wn=wn, btype='lowpass')
    return filtfilt(b, a, sig)
