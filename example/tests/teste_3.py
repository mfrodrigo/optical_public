import numpy as np

from modulator.mach_zehnder_interferometer import MachZehnderInterferometer
from pulse.pulse import Pulse
from pulse.prbs import prbs_7
from output.plotter import Plotter

# dt
from pulse.square_wave import build_square_wave

T = 1000e-12# (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2000  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 3
t = np.arange(0, pulse_number * num_samplesperbit) * dt
t = t - pulse_number * num_samplesperbit * dt / 2
t0 = (np.linspace(-int(pulse_number / 2), int(pulse_number / 2), pulse_number) * T)
t0 = np.repeat(t0, num_samplesperbit)
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=100e-12,
    shift=t0
)
input_seq = np.array([5, 0, 5], ndmin=2).T
squre_wave = build_square_wave(num_samplesperbit=num_samplesperbit, input_seq=input_seq)

Plotter.plot_pulse_input_and_output([[t, abs(pulse.pulse[:, 0]) ** 2, 'pulso original'],
                                    [t, squre_wave[:, 0] / 5000, 'PRBS ']],
                                    graph_title="pulso e onda quadrada",
                                    x_graph='s',
                                    y_graph="")

mach_zehnder = MachZehnderInterferometer(
    pi_voltage_1=5,
    pi_voltage_2=5,
    electro_optical_band=10e9,
    sampling_frequency=1/dt
)
u0 = mach_zehnder.modulate(
    field=pulse.pulse,
    voltage_1=squre_wave,
    voltage_2=5)

Plotter.plot_pulse_input_and_output([[t, np.abs(pulse.pulse[:, 0]) ** 2, "Pulso original"],
                                     [t, np.abs(u0[:, 0]) ** 2, "pulso modulado"],
                                     [t, squre_wave[:, 0] / 5000, 'PRBS ']],
                                    graph_title='Pulso original e modulado',
                                    x_graph='s',
                                    y_graph='w')


