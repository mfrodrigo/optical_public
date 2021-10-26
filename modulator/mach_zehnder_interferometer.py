"""

"""
from numpy import pi, exp
from scipy.signal import butter, lfilter


class MachZehnderInterferometer:
    """

    """

    def __init__(self, pi_voltage_1, pi_voltage_2, electro_optical_band,
                 sampling_frequency):
        self.pi_voltage_1 = pi_voltage_1
        self.pi_voltage_2 = pi_voltage_2
        self.electro_optical_band = electro_optical_band
        self.sampling_frequency = sampling_frequency

    def electro_optical_response(self, input_signal):
        wn = self.electro_optical_band/self.sampling_frequency
        b, a = butter(8, wn, btype='low', analog=False)
        output_signal = lfilter(b, a, input_signal)
        return output_signal

    def modulate(self, field, voltage_1, voltage_2):
        phase_shift_1 = voltage_1 * pi / self.pi_voltage_1
        phase_shift_2 = voltage_2 * pi / self.pi_voltage_2
        output_field = field * 0.5 * (exp(1j * phase_shift_1) + exp(1j * phase_shift_2))
        output_field = self.electro_optical_response(output_field)
        return output_field
