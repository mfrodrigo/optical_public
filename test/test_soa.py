"""
This module test Semiconductor Optical Amplifier (SOA).
"""
import numpy as np
import pytest
from optical_amplifier.soa import SemiconductorOpticalAmplifier


class TestSOA:
    """

    """

    lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
    lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)
    soa = SemiconductorOpticalAmplifier(
        Pin_dbm=0,
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
            self.soa.gain_coefficient(carrier_density=carrier_density)
        assert pytest.approx([material_gain_coefficient[-1],
                              additive_spontaneous_emission_term[-1] / 1e14], abs=1e-1) \
               == [-4.13213e+05, 1.2716]
