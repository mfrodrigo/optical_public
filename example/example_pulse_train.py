import numpy as np
from pulse.pulse import Pulse
from pulse.prbs import prbs_7
from output.plotter import Plotter

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 9  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 128*7
t = (np.array(range(1, pulse_number * num_samplesperbit + 1)) - (pulse_number * num_samplesperbit + 1) / 2) * dt
delta_1 = np.arange(-int((pulse_number+1)/2), 0)*num_samplesperbit+num_samplesperbit/2
delta_2 = np.arange(0, int(pulse_number/2))*num_samplesperbit+num_samplesperbit/2
t0 = np.concatenate((delta_1, delta_2))
t0 = np.repeat(t0, num_samplesperbit) * dt
t = t[0:512*127*7]
t0 = t0[0:512*127*7]
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=100,
    shift=t0
)

array_prbs_7 = np.repeat(prbs_7(), num_samplesperbit, axis=0)/1e3
Plotter.plot_pulse_input_and_output(t, abs(pulse.pulse) ** 2, array_prbs_7, "a")
