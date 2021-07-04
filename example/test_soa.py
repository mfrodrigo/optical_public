"""
TEST SOA
"""
import numpy as np
from matplotlib import pyplot as plt

from optical_amplifier.soa import SemiconductorOpticalAmplifier

t = np.zeros(3200)
t[0] = 0.05e-12
t[-1] = 40
u = np.zeros(3200)
soa = SemiconductorOpticalAmplifier(g0=1096.6/2,
                                    L=6.38336677e-3,
                                    alpha=5)

Pout, pin, gain, phi, tf = soa.amplify_pulse(u, t, 15e-15,
                             FWHM=5e-12)

fig = plt.figure()
plt.plot(tf/1e-12, abs(pin) ** 2, label='input')
plt.plot(tf/1e-12, abs(Pout) ** 2, label='output')
plt.title('Gaussian Pulse ')
plt.ylabel('Power(W)')
plt.xlabel('Time(ps)')
plt.legend()
plt.grid(True)
plt.show()

fig = plt.figure()
plt.plot(tf/1e-12, phi[0:len(t)], label='phase')
plt.title('Gaussian Pulse ')
plt.ylabel('phase rad')
plt.xlabel('Time(ps)')
plt.legend()
plt.grid(True)
plt.show()

fig = plt.figure()
plt.plot(tf/1e-12, gain[0:len(t)], label='output')
plt.title('Gaussian Pulse ')
plt.ylabel('Gain ')
plt.xlabel('Time(ps)')
plt.legend()
plt.grid(True)
plt.show()