"""

"""
import numpy as np
from electrical_components.apply_ltspice import apply_ltspice
from electrical_components.comparator import comparator
from filters.filter_design import bessel_filter
from output.plotter import Plotter

T = 1000e-12
num_samplesperbit = 2000  # should be 2^n
dt = T / num_samplesperbit

t = np.load("/home/rodrigo/Documentos/oitavoperiodo/optical_public/example/tests/time.npy")
i = np.load("/home/rodrigo/Documentos/oitavoperiodo/optical_public/example/tests/corrente-photodetector.npy")
#
t = t - t.min()
configuration_1 = {
    "R": 100
}

time = t
signal_a = i
dummy, signal = apply_ltspice(
    "op_ampt.asc",
    time, signal_a,
    params=configuration_1
)

Plotter.plot_pulse_input_and_output([[time, signal_a, "A"],
                                     [time, signal, "V"]],
                                    graph_title="Resposta do Amaplificador",
                                    x_graph='t',
                                    y_graph="")

sig_filter = bessel_filter(
    1, 10e9, fs=1 / dt, sig=signal
)

Plotter.plot_pulse_input_and_output([[time, signal, "V"],
                                     [time, sig_filter, "V"]],
                                    graph_title="Resposta do filtro Bessel",
                                    x_graph='t',
                                    y_graph="")

output = comparator(-sig_filter, num_samplesperbit, 0.003)

Plotter.plot_pulse_input_and_output([[time, -sig_filter, "V"],
                                     [time, output, "V"]],
                                    graph_title="Resposta do comparador",
                                    x_graph='t',
                                    y_graph="")