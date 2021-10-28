import math
import numpy as np
from pulse.pulse import Pulse
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft
from filter_design import transfer_function

# dt
T = 1000e12
num_samplesperbit = 2 ** 9  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 1
t = (np.array(range(1, pulse_number * num_samplesperbit + 1)) - (pulse_number * num_samplesperbit + 1) / 2) * dt
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=100e12,
)
pulse = pulse.pulse[:, 0]
fig = plt.figure()
plt.plot(t, np.abs(pulse)**2)
plt.xlabel('tempo (s)')
plt.ylabel('gaussian pulse')
plt.grid(True)
fig.savefig("pulso_original.png")

sig_f = fft(pulse)
f = fftfreq(num_samplesperbit, dt)
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

sig_f = fft(sig_filter)
f = fftfreq(num_samplesperbit, 1/dt)
fig = plt.figure()
plt.plot(f, 2.0/num_samplesperbit * np.abs(sig_f))
plt.grid()
plt.xlabel('frequência (Hz)')
plt.ylabel('Gaussian Pulse')
plt.grid(True)
fig.savefig("fft_pulso_filtrado.png")