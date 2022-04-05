"""

"""
import pytest
import numpy as np
from modulator.mach_zehnder_interferometer import MachZehnderInterferometer


class TestMachZehnderInterferometer:

    @pytest.fixture(scope='class')
    def mach_zehnder1(self):
        mach_zehnder = MachZehnderInterferometer(
            pi_voltage_1=5,
            pi_voltage_2=5,
            electro_optical_band=40e9,
            sampling_frequency=80e9
        )
        return mach_zehnder

    @staticmethod
    def test_modulate_destructive(mach_zehnder1):
        input_field = np.ones(100)
        output_field = mach_zehnder1.modulate(
            field=input_field,
            voltage_1=5,
            voltage_2=0)

        np.testing.assert_almost_equal(output_field, np.zeros(100))
