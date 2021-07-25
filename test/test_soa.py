"""
this module test Semiconductor Optical Amplifier (SOA)
"""
import numpy as np
import pytest
from optical_amplifier.soa import SemiconductorOpticalAmplifier


class TestSOA:
    """

    """

    lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
    lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)
    soa = SemiconductorOpticalAmplifier(number_spatial_divisions=100,
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
                              (carrier_density, 'valence', -7.7117e-21)])
    def test_approximation_to_quasi_fermi_level(self, carrier_density, band, expected):
        """
        tests the function that calculates the approximation_to_quasi_fermi_level (nilsson).
        Args:
            carrier_density: (ndarray) Initial guess for carrier density.
            band: (string) band specification (conduction or valence).
            expected: (float) expected value
        """
        answer = self.soa.approximation_to_quasi_fermi_level(carrier_density=carrier_density,
                                                             band=band)
        assert pytest.approx(answer[0], abs=1e-24) == expected

    @pytest.mark.parametrize("carrier_density", [carrier_density])
    def test_gain_coefficient(self, carrier_density):
        material_gain_coefficient, additive_spontaneous_emission_term = \
            self.soa.gain_coefficient(carrier_density=carrier_density)
        assert 1 == 1
