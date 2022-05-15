import numpy as np

from pulse.prbs import prbs_7
from pulse.pulse import Pulse
import matplotlib.pyplot as plt
from filters.filter_design import butter_filter
from output.plotter import Plotter

# dt
T = 1000e-12  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2001  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 127*7
t, t0 = Pulse.generate_time(
    pulse_number, num_samplesperbit,
    dt, T
)
array_prbs_7 = np.repeat(prbs_7(), num_samplesperbit, axis=0) * 5

Plotter.plot_pulse_input_and_output([[t, array_prbs_7[8004:12006, 0] / 5000, 'PRBS ']],
                                    graph_title="Onda quadrada",
                                    x_graph='s',
                                    y_graph="")

sig_filter = butter_filter(
    1, 5e9, fs=1 / dt, sig=array_prbs_7[:, 0]
)

Plotter.plot_pulse_input_and_output([[t, sig_filter[8004:12006] / 5000, 'PRBS ']],
                                    graph_title="Onda filtrada",
                                    x_graph='s',
                                    y_graph="")