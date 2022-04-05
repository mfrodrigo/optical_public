"""
This module has PIN photodetector implementation.
"""
import numpy as np


class PinPhotodiode:
    """
    The PIN photodiode consists in a junction of p-n materials,
    separated by an intrinsic material.
    """
    kb = 1.380649e-23  # J.K-1
    T = 300  # K
    q = 1.602176634e-19  # C

    def __init__(self, quantum_efficiency, wavelength,
                 area, bandwidth, dark_current,
                 load_resistance, noise_figure):
        self.quantum_efficiency = quantum_efficiency
        self.wavelength = wavelength
        self.responsivity = self.calc_responsivity()
        self.area = area
        self.bandwidth = bandwidth
        self.dark_current = dark_current
        self.load_resistance = load_resistance
        self.noise_figure = noise_figure

    def calc_responsivity(self):
        return self.quantum_efficiency * self.wavelength / 1.24

    def calc_shot_noise(self, electric_current):
        return np.sqrt(
            2*self.q*(electric_current+self.dark_current)*self.bandwidth
        )

    def calc_circuit_noise(self):
        return np.sqrt(4 * self.kb * self.T * self.noise_figure * self.bandwidth / self.load_resistance)

    def calc_electric_current(self, electric_field):
        incident_power = np.abs(electric_field) ** 2
