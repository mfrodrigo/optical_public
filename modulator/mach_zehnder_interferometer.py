"""

"""
from numpy import pi, exp, array, linspace, fft
from filters.filter_design import butter_filter


class MachZehnderInterferometer:
    """

    """

    def __init__(self, pi_voltage_1, pi_voltage_2,
                 sampling_frequency,
                 N, cutoff):
        self.pi_voltage_1 = pi_voltage_1
        self.pi_voltage_2 = pi_voltage_2
        self.sampling_frequency = sampling_frequency
        self.N = N
        self.cutoff = cutoff

    def electro_optical_response(self, input_signal, ):
        output_signal = input_signal[:, 0]
        output_signal = butter_filter(self.N, self.cutoff,
                                      self.sampling_frequency,
                                      output_signal)
        return array(output_signal, ndmin=2).T

    def modulate(self, field, voltage_1, voltage_2):
        voltage_1 = self.electro_optical_response(voltage_1)
        phase_shift_1 = voltage_1 * pi / self.pi_voltage_1
        phase_shift_2 = voltage_2 * pi / self.pi_voltage_2
        output_field = field * 0.5 * (exp(1j * phase_shift_1) + exp(1j * phase_shift_2))
        return output_field, voltage_1
