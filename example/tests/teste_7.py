import numpy as np
from channel.channel import Channel
from modulator.mach_zehnder_interferometer import MachZehnderInterferometer
from pulse.pulse import Pulse
from pulse.prbs import prbs_7
from output.plotter import Plotter

# dt
from pulse.square_wave import build_square_wave
from photodetector.pin_photodiode import PinPhotodiode

T = 1000e-12
num_samplesperbit = 2000  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 40

t, t0 = Pulse.generate_time(
    pulse_number, num_samplesperbit,
    dt, T
)
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=40,
    FWHM=100e-12,
    shift=t0
)
input_seq = np.array([5, 0, 5], ndmin=2).T
squre_wave = build_square_wave(num_samplesperbit=num_samplesperbit,
                               size=pulse_number) * 5

Plotter.plot_pulse_input_and_output([[t, abs(pulse.pulse[:, 0]) ** 2, 'pulso original'],
                                     [t, squre_wave[:, 0] / 5000, 'PRBS ']],
                                    graph_title="pulso e onda quadrada",
                                    x_graph='s',
                                    y_graph="")

mach_zehnder = MachZehnderInterferometer(
    pi_voltage_1=5,
    pi_voltage_2=5,
    sampling_frequency=1 / dt,
    N=1,
    cutoff=10e9
)
u0, electro_optical_response = mach_zehnder.modulate(
    field=pulse.pulse,
    voltage_1=squre_wave,
    voltage_2=5)

Plotter.plot_pulse_input_and_output([[t, np.abs(pulse.pulse[:, 0]) ** 2, "Pulso original"],
                                     [t, np.abs(u0[:, 0]) ** 2, "pulso modulado"],
                                     [t, squre_wave[:, 0] / 5000, 'PRBS '],
                                     [t, electro_optical_response[:, 0] / 5000, 'Resposta eletro óptica']],
                                    graph_title='Pulso original e modulado',
                                    x_graph='s',
                                    y_graph='w')

dz = 0.5  # distance stepsize (km)

wavelength = 1550  # nm
speed_of_light = 299792.458  # nm/ps

# Fiber 1
D = 17  # [ps/nm.km]
beta2 = -(D * wavelength ** 2) / (np.pi * speed_of_light)  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial
gamma = 0.01
alpha = 0.2 / 4.343

# DCE Fiber
D_DCE = -100  # [ps/nm.km]
beta2_DCE = -(D_DCE * (wavelength ** 2)) / (np.pi * speed_of_light)  # beta2 (ps^2/km)
betap_DCE = np.transpose(np.array([0, 0, beta2_DCE]).reshape(1, 3))  # dispersion polynomial
gamma_DCE = 0.03
alpha_DCE = 0.4 / 4.343

nz_step = [10]

for nz in nz_step:
    # output
    u1 = Channel.ssprop(u0, dt * 1e12, dz, nz, alpha, betap, gamma)

    title_graph_1 = 'Plot canal: ' + str(nz * dz) + 'Km, alpha = ' \
                    + str(alpha) + '_beta_2_' + str(beta2) + '_gamma_' \
                    + str(gamma) + '.png'

    Plotter.plot_pulse_input_and_output([[t, abs(u0[:, 0]) ** 2, "Input da fibra"],
                                         [t, abs(u1[:, 0]) ** 2, "Output da fibra"]],
                                        graph_title=title_graph_1,
                                        x_graph='s',
                                        y_graph="")

    nz_DCE = -int(nz * D / D_DCE)
    u2 = Channel.ssprop(u1, dt * 1e12, dz, nz_DCE, alpha_DCE, betap_DCE, gamma_DCE)

    Plotter.plot_pulse_input_and_output([[t, np.abs(u1[:, 0]) ** 2, "Input DCF"],
                                         [t, np.abs(u2[:, 0]) ** 2, "Output DCF"]],
                                        graph_title="saida dcf",
                                        x_graph='s',
                                        y_graph="")

pin = PinPhotodiode(
    quantum_efficiency=0.8,
    wavelength=1.550,  # 1530 nm
    bandwidth=5e9,
    dark_current=1e-9,  # nA
    load_resistance=1000,
    noise_figure_db=3
)

i = pin.calc_electric_current(u2[:, 0])
Plotter.plot_pulse_input_and_output([[t, i, "A"]],
                                    graph_title="Resposta do fotodiodo Pin",
                                    x_graph='t',
                                    y_graph="")
