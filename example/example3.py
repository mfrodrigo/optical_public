"""

"""
import math
import numpy as np
from channel.channel import Channel
from pulse.pulse import Pulse
from pulse.half_power import return_half_power
from output.plotter import Plotter
from optical_amplifier.soa import SemiconductorOpticalAmplifier
from output.tables import Tables

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 9  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=100
)

dz = 0.5  # distance stepsize (km)

wavelength = 1550  # nm
speed_of_light = 299792.458  # nm/ps

# Fiber 1
D = 17  # [ps/nm.km]
beta2 = -(D * wavelength ** 2) / (math.pi * speed_of_light)  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial
gamma = 0.01
alpha = 0.2 / 4.343

# DCE Fiber
D_DCE = -100  # [ps/nm.km]
beta2_DCE = -(D_DCE * (wavelength ** 2)) / (math.pi * speed_of_light)  # beta2 (ps^2/km)
betap_DCE = np.transpose(np.array([0, 0, beta2_DCE]).reshape(1, 3))  # dispersion polynomial
gamma_DCE = 0.03
alpha_DCE = 0.4 / 4.343

nz_step = [10]

lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)

for nz in nz_step:
    # output
    u1 = Channel.ssprop(pulse.pulse, dt, dz, nz, alpha, betap, gamma)
    pulse.original_pulse = Channel.ssprop(pulse.original_pulse, dt, dz, nz, alpha, betap, gamma)
    title_graph_1 = 'Plot canal: ' + str(nz * dz) + 'Km, alpha = ' \
                    + str(alpha) + '_beta_2_' + str(beta2) + '_gamma_' \
                    + str(gamma) + '.png'

    Plotter.plot_pulse_input_and_output([[t, abs(pulse.pulse[:, 0]) ** 2, "Input da fibra"],
                                         [t, abs(u1[:, 0]) ** 2, "Output da fibra"]],
                                        graph_title=title_graph_1,
                                        x_graph='s',
                                        y_graph="")

    nz_DCE = -int(nz * D / D_DCE)
    u2 = Channel.ssprop(u1, dt, dz, nz_DCE, alpha_DCE, betap_DCE, gamma_DCE)
    pulse.original_pulse = Channel.ssprop(pulse.original_pulse, dt, dz, nz_DCE, alpha_DCE, betap_DCE, gamma_DCE)

    Plotter.plot_pulse_input_and_output([[t, np.abs(u1[:, 0]) ** 2, "Input DCF"],
                                         [t, np.abs(u2[:, 0]) ** 2, "Output DCF"]],
                                        graph_title="saida dcf",
                                        x_graph='s',
                                        y_graph="")

