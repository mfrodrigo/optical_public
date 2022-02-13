import numpy as np
from pulse.pulse import Pulse
import matplotlib.pyplot as plt
from filters.filter_design import butter_filter
from output.plotter import Plotter

# dt
T = 1000e-12
num_samplesperbit = 2001  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 1
t = np.arange(0, pulse_number * num_samplesperbit) * dt
t = t - pulse_number * num_samplesperbit * dt / 2
pulse = Pulse(
    power=1,
    time=t,
    SNR_dB=100,
    FWHM=400e-12,
)
pulse = pulse.pulse[:, 0]

Plotter.plot_pulse_input_and_output([[t, np.abs(pulse) ** 2, " pulso original"]],
                                    graph_title="Pulso original no domínio do tempo",
                                    x_graph="s",
                                    y_graph='Watt'
                                    )

sig_f = np.fft.fftshift(np.fft.fft(pulse)) / len(pulse)
f = np.linspace(-0.5, 0.5, len(pulse)) * 1 / dt

Plotter.plot_pulse_input_and_output([[f, np.abs(sig_f) ** 2, " FFT pulso original"]],
                                    graph_title="FFT pulso original",
                                    x_graph="Hz",
                                    y_graph='Watt'
                                    )

sig_filter = np.fft.ifftshift(np.fft.ifft(sig_f * butter_filter(
    4, 1e9, fs=1 / dt, sig=f
))) * len(sig_f)
fig = plt.figure()
plt.plot(t, np.abs(sig_filter) ** 2,
         label='Resposta ao filtro')
plt.xlabel('tempo (s)')
plt.ylabel('gaussian pulse')
plt.legend()
plt.grid(True)
fig.savefig("pulso_filtrado.png")

sig_f_output = np.fft.fft(np.fft.fftshift(sig_filter)) / len(sig_filter)
fig = plt.figure()
plt.plot(f, np.abs(sig_f_output))
plt.grid()
plt.xlabel('frequência (Hz)')
plt.ylabel('Gaussian Pulse')
plt.grid(True)
fig.savefig("fft_pulso_filtrado.png")
#
#
Plotter.plot_pulse_input_and_output([[t, np.abs(pulse) ** 2, "input"],
                                     [t, np.abs(sig_filter) ** 2, "output"]],
                                    graph_title="Comparação dos pulsos de entrada e saída",
                                    x_graph='s',
                                    y_graph='Watt')
Plotter.plot_pulse_input_and_output([[f, np.abs(sig_f) ** 2, 'input'],
                                     [f, np.abs(sig_f_output) ** 2, 'output']],
                                    graph_title="FFT dos pulsos de entrada e saída",
                                    x_graph='Hz',
                                    y_graph='Watt')

# Plotter.plot_pulse_input_and_output([[f, 20 * np.log10(np.abs(transfer_function(f))), 'filtro']],
#                                     graph_title="Função de transferência",
#                                     x_graph='Hz',
#                                     y_graph='dB'
#                                     )
