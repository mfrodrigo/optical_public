"""
Test Soa
"""
import numpy as np
from optical_amplifier.soa import SemiconductorOpticalAmplifier

lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)
soa = SemiconductorOpticalAmplifier(
        Pin_dbm=np.arange(-40, 15, 5),
        wavelength_s=1550,
        number_spatial_divisions=100,
        number_spectrum_slices=100,
        wavelength_0=lambda0,
        wavelength_1=lambda1)

Pout_dBm, Gain, noise_figure = soa.run_simulation_soa()
print(Pout_dBm)
print(Gain)
print(noise_figure)