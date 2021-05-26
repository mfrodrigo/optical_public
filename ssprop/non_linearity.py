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

# dt
T = 1000  #  (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 11  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

dz = 0.5  # distance stepsize (km)
nz = 40

beta2 = -43  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial

t0 = 0
C = 0
m = 1
P0 = 0.1

##
n2 = 2.6*(10**-20)
gamma = 2*math.pi*(10**24)*n2/(1550*76)

gamma = 0.0
alpha = 0.0
# u0
FWHM = 100
u0 = np.zeros(shape=(len(t), 1), dtype=complex)
u0[:, 0] = math.sqrt(P0) * 2 ** (-((1 + 1j * C) / 2) * (2 * (t - t0) / FWHM) ** (2 * m))
# Output

u1 = np.zeros(shape=(len(t), nz), dtype=complex)
u1[:, 0] = u0[:, 0]

for i in range(1, nz):
    aux = u1[:, i - 1]
    aux = aux.reshape(len(t), 1)
    u1[:, i] = ssprop(aux, dt, dz, 1, alpha, betap, gamma)

print('###########################################################')
u_output = abs(u1[:, nz-1]) ** 2
max_peak = np.amax(u_output)
for i in range(u_output.shape[0]):
    if u_output[i] >= max_peak/2:
        aux_1 = i
        break
for j in range(aux_1, u_output.shape[0]):
    if u_output[j] <= max_peak/2:
        aux_2 = j
        break

delta_t = t[aux_2]-t[aux_1]
print(t[aux_1], t[aux_2], delta_t, max_peak, u_output[aux_1], u_output[aux_2])
print('###########################################################')
plt.figure()

plt.plot(t, abs(u1[:, 0]) ** 2, label=f'$u_{0}$')
plt.plot(t, abs(u1[:, nz-1]) ** 2, label=f'$u_{nz}$')
plt.title('Gaussian Pulse ')
plt.xlabel(r'$(t-\beta_1z)/T_0$')
plt.ylabel('|u1(z,t)|^2/P_0')
plt.legend()
plt.grid(True)
plt.show()

