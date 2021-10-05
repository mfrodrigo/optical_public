import math
import numpy as np
from pulse.pulse import Pulse
from pulse.prbs import prbs_7
from channel.channel import Channel
from output.plotter import Plotter
from modulator.mach_zehnder_interferometer import MachZehnderInterferometer

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 9  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 128 * 7
t = (np.array(range(1, pulse_number * num_samplesperbit + 1)) - (pulse_number * num_samplesperbit + 1) / 2) * dt
delta_1 = np.arange(-int((pulse_number + 1) / 2), 0) * num_samplesperbit + num_samplesperbit / 2
delta_2 = np.arange(0, int(pulse_number / 2)) * num_samplesperbit + num_samplesperbit / 2
t0 = np.concatenate((delta_1, delta_2))
t0 = np.repeat(t0, num_samplesperbit) * dt
t = t[0:512 * 127 * 7]
t0 = t0[0:512 * 127 * 7]
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=100,
    shift=t0
)

array_prbs_7 = np.repeat(prbs_7(), num_samplesperbit, axis=0) * 5
Plotter.plot_pulse_input_and_output(t, abs(pulse.pulse) ** 2, array_prbs_7, "a")

mach_zehnder = MachZehnderInterferometer(
    pi_voltage_1=5,
    pi_voltage_2=5,
    electro_optical_band=40e9,
    sampling_frequency=80e9
)
u0 = mach_zehnder.modulate(
    field=pulse.pulse,
    voltage_1=array_prbs_7,
    voltage_2=5)

Plotter.plot_pulse_input_and_output(t, abs(pulse.pulse) ** 2, abs(u0) ** 2, "a")
Plotter.plot_pulse_input_and_output(t, abs(u0) ** 2, array_prbs_7 / 5000, "a")

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

nz_step = [170]

lambda0 = 1300  # start wavelength for gain coefficient and ASE spectrum (nm)
lambda1 = 1650  # end wavelength for gain coefficient and ASE spectrum (nm)

for nz in nz_step:
    # output
    u1 = Channel.ssprop(u0, dt, dz, nz, alpha, betap, gamma)
    title_graph_1 = 'Plot canal: ' + str(nz * dz) + 'Km, alpha = ' \
                    + str(alpha) + '_beta_2_' + str(beta2) + '_gamma_' \
                    + str(gamma) + '.png'

    Plotter.plot_pulse_input_and_output(t, abs(u0) ** 2, abs(u1) ** 2, title_graph_1)

    nz_DCE = -int(nz * D / D_DCE)
    u2 = Channel.ssprop(u1, dt, dz, nz_DCE, alpha_DCE, betap_DCE, gamma_DCE)
    Plotter.plot_pulse_input_and_output(t, abs(u1) ** 2, abs(u2) ** 2, "Fibra DCF")

    title_graph_1 = 'Input_output.png'

    Plotter.plot_pulse_input_and_output(t, abs(u0) ** 2, abs(u2) ** 2, title_graph_1)
