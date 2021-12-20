import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft
from filters.filter_design import transfer_function
from pulse.pulse import Pulse
from pulse.prbs import prbs_7
from output.plotter import Plotter

# dt
T = 1000e-12  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 9  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 127*7
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

array_prbs_7 = np.repeat(prbs_7(), num_samplesperbit, axis=0) * 5
Plotter.plot_pulse_input_and_output(t, abs(pulse.pulse)[:, 0] ** 2, array_prbs_7[:, 0] / 5000, "a")

pulse = pulse.pulse[:, 0]
fig = plt.figure()
plt.plot(t, np.abs(pulse)**2)
plt.xlabel('tempo (s)')
plt.ylabel('gaussian pulse')
plt.grid(True)
fig.savefig("pulso_original.png")

sig_f = fft(pulse)
f = fftfreq(num_samplesperbit*pulse_number, dt)
fig = plt.figure()
plt.plot(f, 2.0/num_samplesperbit * np.abs(sig_f))
plt.grid()
plt.xlabel('frequência (Hz)')
plt.ylabel('Gaussian Pulse')
plt.grid(True)
fig.savefig("fft_gausssian_pulse.png")

sig_filter = ifft(sig_f*transfer_function(f))
fig = plt.figure()
plt.plot(t, np.abs(sig_filter)**2, label='Resposta ao filtro')
plt.xlabel('tempo (s)')
plt.ylabel('gaussian pulse')
plt.legend()
plt.grid(True)
fig.savefig("pulso_filtrado.png")

sig_f_output = fft(sig_filter)
f = fftfreq(num_samplesperbit*pulse_number, dt)
fig = plt.figure()
plt.plot(f, 2.0/num_samplesperbit * np.abs(sig_f))
plt.grid()
plt.xlabel('frequência (Hz)')
plt.ylabel('Gaussian Pulse')
plt.grid(True)
fig.savefig("fft_pulso_filtrado.png")

Plotter.plot_pulse_input_and_output(t, abs(pulse) ** 2, abs(sig_filter) ** 2, "Pulsos")
Plotter.plot_pulse_input_and_output(f, abs(sig_f), abs(sig_f_output), "Pulsos")
