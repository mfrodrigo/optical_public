"""
Pulse Class
"""
import numpy as np
from math import sqrt


class Pulse:

    def __init__(self, power, time, SNR_dB,
                 FWHM, complex_factor=0,
                 shift=0, exponent=1,
                 type_pulse="Gaussian"):
        """

        Args:
            power:
            time:
            SNR_dB:
            FWHM:
            complex_factor:
            shift:
            exponent:
            type_pulse:
        """
        self.P0 = power
        self.t = time
        self.SNR = 10 ** (SNR_dB / 10)
        self.FWHM = FWHM
        self.C = complex_factor
        self.t0 = shift
        self.m = exponent
        self.type_pulse = type_pulse
        self.pulse = None
        if type_pulse == "Gaussian":
            self._gaussian_pulse()
        self.original_pulse = self.pulse.copy()
        self._add_noise()

    def _gaussian_pulse(self):
        """
        """
        self.pulse = np.zeros(shape=(len(self.t), 1), dtype=complex)
        self.pulse[:, 0] = sqrt(self.P0) * 2 ** (
                    -((1 + 1j * self.C) / 2) * (2 * (self.t - self.t0) / self.FWHM) ** (2 * self.m))

    def _add_noise(self):
        """
        Returns:

        """
        self.pulse[:, 0] = np.sqrt(abs(self.pulse[:, 0])**2 + self.P0 / self.SNR * np.random.randn(self.pulse.shape[0]),
                                   dtype="complex")
