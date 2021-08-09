import math
import numpy as np
from channel.channel import Channel
from pulse.pulse import Pulse
from pulse.half_power import return_half_power
from output.plotter import Plotter
from optical_amplifier.soa import SemiconductorOpticalAmplifier
from output.tables import Tables


T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 11  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

pulse_1 = Pulse(
    power=100e-6,
    time=t,
    SNR_dB=0,
    FWHM=100
)

pulse_2 = Pulse(
    power=100e-6,
    time=t,
    SNR_dB=2000,
    FWHM=100
)

Plotter.plot_pulse_input_and_output(t, abs(pulse_1.pulse) ** 2, abs(pulse_2.pulse) ** 2, "pulso ruidoso")
Plotter.plot_pulse_input_and_output(t, abs(pulse_1.pulse) ** 2-abs(pulse_2.pulse) ** 2, abs(pulse_1.pulse) ** 2-abs(pulse_2.pulse) ** 2, "pulso ruidoso 2")