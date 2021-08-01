"""
This module test Semiconductor Optical Amplifier (SOA).
"""
from math import sqrt
import numpy as np
import pytest
from optical_amplifier.soa import SemiconductorOpticalAmplifier


class TestSOA:
    """

    """

    lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
    lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)
    soa = SemiconductorOpticalAmplifier(
        Pin_dbm=np.arange(-40, 15, 5),
        wavelength_s=1550,
        number_spatial_divisions=100,
        number_spectrum_slices=100,
        wavelength_0=lambda0,
        wavelength_1=lambda1)
    carrier_density = np.ones(100) * 1.2e24

    @pytest.mark.parametrize('carrier_density', [carrier_density])
    def test_energy_gap(self, carrier_density):
        """
        tests the function that calculates the energy gap.
        Args:
            carrier_density: (ndarray) Initial guess for carrier density.
        """
        expected = 122.83044e-19
        assert pytest.approx(np.sum(self.soa.energy_gap(carrier_density=carrier_density)),
                             abs=1e-20) == expected

    @pytest.mark.parametrize("carrier_density, band, expected",
                             [(carrier_density, 'conduction', 13.5808e-21),
                              (carrier_density, 'valence', 7.7117e-21)])
    def test_approximation_to_quasi_fermi_level(self, carrier_density, band, expected):
        """
        Tests the function that calculates the approximation_to_quasi_fermi_level (nilsson).
        Args:
            carrier_density: (ndarray) Initial guess for carrier density.
            band: (string) band specification (conduction or valence).
            expected: (float) expected value
        """
        answer = self.soa.approximation_to_quasi_fermi_level(carrier_density=carrier_density,
                                                             band=band)
        assert pytest.approx(answer[0], abs=1e-24) == expected

    def test_gain_coefficient(self):
        """Tests the function that calculates the material gain coefficient and
         additive spontaneous emission term"""
        carrier_density = np.ones(101) * 1.2e24
        material_gain_coefficient, additive_spontaneous_emission_term = \
            self.soa.gain_coefficient(carrier_density=carrier_density,
                                      energy=self.soa.energy)
        assert pytest.approx([material_gain_coefficient[-1],
                              additive_spontaneous_emission_term[-1] / 1e14], abs=1e-1) \
               == [-4.13213e+05, 1.2716]

    @pytest.mark.parametrize('carrier_density', [carrier_density])
    def test_calc_alpha(self, carrier_density):
        """Tests calc_alpha function."""
        alpha_s, alpha = self.soa.calc_alpha(carrier_density=carrier_density)
        assert [alpha_s[99], alpha[99, 100]] == [10250, 10250]

    @pytest.mark.parametrize('carrier_density', [carrier_density])
    def test_solve_travelling_wave_equations_ASE(self, carrier_density):
        alpha_s, _ = self.soa.calc_alpha(carrier_density=carrier_density)
        number_division = self.soa.number_spatial_divisions
        forward_signal_amplitude = np.zeros(number_division + 1, dtype=np.complex128)
        backward_signal_amplitude = np.zeros(number_division + 1, dtype=np.complex128)
        forward_signal_amplitude[0] = (1 - self.soa.r1) * \
                                      sqrt(self.soa.eta_in) * 8.8277e5 + \
                                      self.soa.r1 * backward_signal_amplitude[0]
        backward_signal_amplitude[-1] = self.soa.r2 * forward_signal_amplitude[-1]
        material_gain_coefficient_signal, _ = self.soa.gain_coefficient(
            carrier_density=carrier_density,
            energy=self.soa.energy_signal)
        material_gain_coefficient_signal = material_gain_coefficient_signal[0:self.soa.number_spatial_divisions]
        forward_signal_amplitude, backward_signal_amplitude = self.soa.solve_travelling_wave_equations_signal(
            forward_signal_amplitude, backward_signal_amplitude,
            material_gain_coefficient_signal, alpha_s
        )

        assert pytest.approx([forward_signal_amplitude[99], backward_signal_amplitude[10]], abs=1)\
               == [-1.3834e+04 - 2.2584e+04j, 0]

    def test_run_simulation_soa(self):
        self.soa.run_simulation_soa()
