# test ssprop

# INPUT
#
# u0 - starting field amplitude (vector)
# dt - time step
# dz - propagation stepsize
# nz - number of steps to take, ie, ztotal = dz*nz
# alpha - power loss coefficient, ie, P=P0*exp(-alpha*z)
# betap - dispersion polynomial coefs, [beta_0 ... beta_m]
# gamma - nonlinearity coefficient
# maxiter - max number of iterations (default = 4)
# tol - convergence tolerance (default = 1e-5)
#
# OUTPUT
#
# u1 - field at the output

# Libraries

import math
import numpy as np
from scipy import signal
from ssprop import ssprop
import matplotlib.pyplot as plt
from scipy.stats import hypsecant

# dt
T = 32  # pulse width FWHM (ps)
num_samplesperbit = 2 ** 10  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

dz = 2  # distance stepsize (km)
nz = 2

beta2 = 1  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial

t0 = 0
C = 0
m = 1
P0 = 1

# u0
FWHM = 2 * math.sqrt(math.log(2))
u0 = np.zeros(shape=(len(t), 1), dtype=complex)
u0[:, 0] = hypsecant.pdf(t - t0)
# Output

u1 = np.zeros(shape=(len(t), nz - 1), dtype=complex)

for i in range(nz - 1):
    aux = u0[:, 0]
    aux = aux.reshape(len(t), 1)
    u1[:, i] = ssprop(aux, dt, dz, 1, 0, betap, 0)

print('###########################################################')
print(u1)
plt.figure()
plt.plot(t, abs(u0) ** 2, label=r'$u_o$')
plt.plot(t, abs(u1) ** 2, label=r'$u_1$')
plt.title('Hyperbolic Secant')
plt.xlabel(r'$(t-\beta_1z)/T_0$')
plt.ylabel('|u1(z,t)|^2/P_0')
plt.legend()
plt.grid(True)
plt.show()
