import math
import numpy as np
from pulse.pulse import Pulse
from scipy import signal
from pulse.prbs import prbs_7
from channel.channel import Channel
from output.plotter import Plotter
import matplotlib.pyplot as plt
from modulator.mach_zehnder_interferometer import MachZehnderInterferometer

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 10  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 1
t = (np.array(range(1, pulse_number * num_samplesperbit + 1)) - (pulse_number * num_samplesperbit + 1) / 2) * dt
delta_1 = np.arange(-int((pulse_number + 1) / 2), 0) * num_samplesperbit + num_samplesperbit / 2
delta_2 = np.arange(0, int(pulse_number / 2)) * num_samplesperbit + num_samplesperbit / 2
t0 = np.concatenate((delta_1, delta_2))
t0 = np.repeat(t0, num_samplesperbit) * dt
t0 = t0[0:512]
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=0.33*T,
    shift=0
)

# array_prbs_7 = np.repeat(prbs_7(1), num_samplesperbit, axis=0) * 5
Plotter.plot_pulse_input_and_output(t, abs(pulse.pulse) ** 2, abs(pulse.pulse) ** 2, "a")

sos = signal.butter(4, 2*np.pi*40e9, 'low', fs=(1/dt)*1e12, output='sos')
filtered = signal.sosfilt(sos, pulse.pulse)
Plotter.plot_pulse_input_and_output(1/t, abs(pulse.pulse) ** 2, abs(np.fft.fft(pulse.pulse, axis=0)), "a")

b, a = signal.butter(4, 40e9/1e12, 'low', analog=False)
w, h = signal.freqs(b, a)
plt.semilogx(w, 20 * np.log10(abs(h)))
plt.title('Butterworth filter frequency response')
plt.xlabel('Frequency [radians / second]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(100, color='green') # cutoff frequency
plt.show()