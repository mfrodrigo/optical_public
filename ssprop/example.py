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
from half_power import return_half_power
import matplotlib.pyplot as plt

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 11  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

dz = 0.5  # distance stepsize (km)
nz = 120

beta2 = -43  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial

t0 = 0
C = 0
m = 1
P0 = 0.1

##
n2 = 2.6 * (10 ** -20)
# gamma = 2*math.pi*(10**24)*n2/(1550*76)

gamma = 0.01
alpha = 0.0
# u0
FWHM = 100
u0 = np.zeros(shape=(len(t), 1), dtype=complex)
u0[:, 0] = math.sqrt(P0) * 2 ** (-((1 + 1j * C) / 2) * (2 * (t - t0) / FWHM) ** (2 * m))

# output
u1 = ssprop(u0, dt, dz, nz, alpha, betap, gamma)

print('###########################################################')
list_values = return_half_power(t, u1)
print(list_values)
print('###########################################################')
fig = plt.figure()
plt.plot(t, abs(u0[:, 0]) ** 2, label='Input')
plt.plot(t, abs(u1[:, 0]) ** 2, label='Output')
plt.title('Gaussian Pulse ')
plt.xlabel(r'$(t-\beta_1z)/T_0$')
plt.ylabel('|u1(z,t)|^2/P_0')
plt.legend()
plt.grid(True)
fig.savefig('Plot canal: ' + str(nz * dz) + 'Km, alpha = ' + str(alpha) + '.png', dpi=fig.dpi)
plt.show()