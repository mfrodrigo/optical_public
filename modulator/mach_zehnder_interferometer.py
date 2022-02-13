"""

"""
from numpy import pi, exp, array, linspace, fft
from filters.filter_design import transfer_function


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
        output_signal = input_signal[:, 0]
        f = linspace(-0.5, 0.5, len(output_signal)) * self.sampling_frequency
        output_signal = fft.fftshift(fft.fft(output_signal)) / len(output_signal)
        output_signal = fft.ifftshift(fft.ifft(output_signal * transfer_function(f))) * len(output_signal)
        return array(output_signal, ndmin=2).T

    def modulate(self, field, voltage_1, voltage_2):
        field = self.electro_optical_response(voltage_1)
        phase_shift_1 = voltage_1 * pi / self.pi_voltage_1
        phase_shift_2 = voltage_2 * pi / self.pi_voltage_2
        output_field = field * 0.5 * (exp(1j * phase_shift_1) + exp(1j * phase_shift_2))
        return output_field
