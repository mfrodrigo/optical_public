"""

"""
import numpy as np
from apply_ltspice import apply_ltspice
from output.plotter import Plotter

t = np.load("/home/rodrigo/Documentos/oitavoperiodo/optical_public/electrical_components/time.npy")
i = np.load("/home/rodrigo/Documentos/oitavoperiodo/optical_public/electrical_components/corrente-photodetector.npy")
#
t = t - t.min()
# configuration_1 = {
#     "R": 100
# }

# our samples shall be 100 ms wide
# time = np.linspace(0, 1, int(1e4), False)  # 1 second
# signal_a = 1e-6*np.sin(2*np.pi*100*time)

##################################################
##        apply filter - configuration 1        ##
##################################################
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