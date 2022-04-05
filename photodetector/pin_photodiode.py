"""
This module has PIN photodetector implementation.
"""
import numpy as np

class PinPhotodiode:
    """
    The PIN photodiode consists in a junction of p-n materials,
    separated by an intrinsic material.
    """

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

    def calc_shot_noise(self):
        pass

    def calc_circuit_noise(self):
        pass

    def calc_electric_current(self, electric_field):
        incident_power = np.abs(electric_field) ** 2