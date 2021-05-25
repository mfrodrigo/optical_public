"""

"""
import math
import numpy as np 

# CONSTANTS
c = 299792458  # speed of light (m/s)

# NUMERICAL PARAMETERS
numbitspersymbol = 1
P0 = 0.003  # peak power (W)
FWHM = 25  # pulse width FWHM (ps)

nt = 2 ** 8  # number of points in FFT
dt = FWHM / nt  # sampling time(ps)

dz = 0.2  # distance stepsize (km)
nz = 500  # step number

wavelength = 1550  # wavelength (nm)
alpha_indB = 0.17  # loss (dB/km)
D = 18.5  # GVD (ps/nm.km)
beta3 = 0.06  # GVD slope (ps^3/km)

# CALCULATED QUANTITIES
alpha_loss = math.log(math.e, 10)*alpha_indB/10
beta2 = -1000*D*wavelength**2/(2*math.pi*c)

betap = np.transpose(np.array([0, 0, beta2, beta3]).reshape(1, 4))  # dispersion polynomial

