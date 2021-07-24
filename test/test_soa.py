"""
this module test Semiconductor Optical Amplifier (SOA)
"""
import numpy as np
import pytest
from optical_amplifier.soa import SemiconductorOpticalAmplifier


class TestSOA:
    """

    """
    soa = SemiconductorOpticalAmplifier()

    carrier_density = np.ones(100) * 1.2e24

    @pytest.mark.parametrize('carrier_density', [carrier_density])
    def test_energy_gap(self, carrier_density):
        """
        tests the function that calculates the energy gap.
        Args:
            carrier_density: (ndarray)
        """
        expected = 122.83044e-19
        assert pytest.approx(np.sum(self.soa.energy_gap(carrier_density=carrier_density)),
                             abs=1e-20) == expected

    @pytest.mark.parametrize("carrier_density, band, expected",
                             [(carrier_density, 'conduction', 13.5808e-21),
                              (carrier_density, 'valence', -7.7117e-21)])
    def test_approximation_to_quasi_fermi_level(self, carrier_density, band, expected):
        answer = self.soa.approximation_to_quasi_fermi_level(carrier_density=carrier_density,
                                                             band=band)
        assert pytest.approx(answer[0], abs=1e-24) == expected
