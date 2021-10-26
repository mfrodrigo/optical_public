import math
import numpy as np
from pulse.pulse import Pulse
from scipy import signal
from pulse.prbs import prbs_7
from channel.channel import Channel
from scipy.fft import fft, fftfreq
import plotly.graph_objects as go
from output.plotter import Plotter
import matplotlib.pyplot as plt
from modulator.mach_zehnder_interferometer import MachZehnderInterferometer
from scipy.signal import butter, lfilter

# dt
T = 6000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 9  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
pulse_number = 1
t = (np.array(range(1, pulse_number * num_samplesperbit + 1)) - (pulse_number * num_samplesperbit + 1) / 2) * dt
delta_1 = np.arange(-int((pulse_number + 1) / 2), 0) * num_samplesperbit + num_samplesperbit / 2
delta_2 = np.arange(0, int(pulse_number / 2)) * num_samplesperbit + num_samplesperbit / 2
t0 = np.concatenate((delta_1, delta_2))
t0 = np.repeat(t0, num_samplesperbit) * dt
t = t[0:512]*1e-12
t0 = t0[0:512]*1e-12
pulse = Pulse(
    power=1e-3,
    time=t,
    SNR_dB=30,
    FWHM=1000*1e-12,
    shift=t0
)

# array_prbs_7 = np.array(np.repeat(prbs_7(), num_samplesperbit, axis=0), dtype=complex)
# pulse.pulse = array_prbs_7*pulse.pulse
# fig = go.Figure()
# fig.add_trace(go.Scatter(
#             x=t,
#             y=np.abs(pulse.pulse/np.max(np.abs(pulse.pulse)))[:, 0],
#         ))
# fig.show()

yf = fft(pulse.pulse)
xf = fftfreq(num_samplesperbit, dt)[:num_samplesperbit//2]
fig = go.Figure()
fig.add_trace(go.Scatter(
            x=xf,
            y=np.abs(yf[:num_samplesperbit//2, 0]),
        ))
fig.show()

sos = signal.butter(4, 2*np.pi*10, 'low', fs=2*np.pi*22, output='sos')
output_signal = signal.sosfilt(sos, pulse.pulse/np.max(np.abs(pulse.pulse)))
output_signal = output_signal/np.max(np.abs(output_signal))
fig = go.Figure()
fig.add_trace(go.Scatter(
            x=t,
            y=np.abs(output_signal[:, 0])**2,
            name='output',
        ))
fig.add_trace(go.Scatter(
            x=t,
            y=(np.abs(pulse.pulse/np.max(np.abs(pulse.pulse)))[:, 0]) **2,
            name='input',
        ))
fig.show()

print("OI")