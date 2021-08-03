# test channel

# INPUT
#
# u0 - starting field amplitude (vector)
# dt - time step
# dz - propagation stepsize
# nz - number of steps to take, ie, ztotal = dz*nz
# alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
# betap - dispersion polynomial coefs, [beta_0 ... beta_m]
# gamma - nonlinearity coefficient
# maxiter - max number of iterations (default = 4)
# tol - convergence tolerance (default = 1e-5)
#
# OUTPUT
#
# u1 - field at the output

# Libraries

import math
import numpy as np
from channel.channel import Channel
from pulse.half_power import return_half_power
from output.plotter import Plotter
from optical_amplifier.soa import SemiconductorOpticalAmplifier
from output.tables import Tables

# dt
T = 500 # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 7  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

# Pulse
FWHM = 100
t0 = 0
C = 0
m = 1
P0 = 0.01
u0 = np.zeros(shape=(len(t), 1), dtype=complex)
u0[:, 0] = math.sqrt(P0) * 2 ** (-((1 + 1j * C) / 2) * (2 * (t - t0) / FWHM) ** (2 * m))
dz = 0.5  # distance stepsize (km)

wavelength = 1550  # nm
speed_of_light = 299792.458  # nm/ps

# Fiber 1
D = 17  # [ps/nm.km]
beta2 = -(D * wavelength ** 2) / (math.pi * speed_of_light)  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial
gamma = 0.0
alpha = 0.2/4.343

# DCE Fiber
D_DCE = -100  # [ps/nm.km]
beta2_DCE = -(D_DCE * (wavelength ** 2)) / (math.pi * speed_of_light)  # beta2 (ps^2/km)
betap_DCE = np.transpose(np.array([0, 0, beta2_DCE]).reshape(1, 3))  # dispersion polynomial
gamma_DCE = 0.03
alpha_DCE = 0.4/4.343

nz_step = [10]

lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)

for nz in nz_step:
    # output
    u1 = Channel.ssprop(u0, dt, dz, nz, alpha, betap, gamma)

    title_graph_1 = 'Plot canal: ' + str(nz * dz) + 'Km, alpha = ' \
                    + str(alpha) + '_beta_2_' + str(beta2) + '_gamma_' \
                    + str(gamma) + '.png'

    Plotter.plot_pulse_input_and_output(t, abs(u0)**2, abs(u1)**2, title_graph_1)

    u3 = 10*np.log10(abs(u1)**2/1e-3).transpose()[0]
    soa = SemiconductorOpticalAmplifier(
        Pin_dbm=u3,
        wavelength_s=wavelength,
        number_spatial_divisions=100,
        number_spectrum_slices=100,
        wavelength_0=lambda0,
        wavelength_1=lambda1)

    Pout_dBm, Gain, noise_figure = soa.run_simulation_soa()

    u4 = (abs(u1.transpose())**2)*Gain
    title_graph_1 = 'Soa.png'

    Plotter.plot_pulse_input_and_output(t, abs(u1)**2, u4[0], title_graph_1)

    title_graph_1 = 'Input_output.png'

    Plotter.plot_pulse_input_and_output(t,  abs(u0)**2, u4[0], title_graph_1)