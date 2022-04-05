# test channel

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
import pandas as pd
from scipy import signal
from channel.channel import Channel
from output.plotter import Plotter
from pulse.half_power import return_half_power
import matplotlib.pyplot as plt

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 11  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

dz = 0.5  # distance stepsize (km)
nz = 160

beta2 = -43  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial

t0 = 0
C = 0
m = 1
P0 = 1

##
n2 = 2.6 * (10 ** -20)
# gamma = 2*math.pi*(10**24)*n2/(1550*76)

gamma = 0.01
alpha = 0.
# u0
FWHM = 100
u0 = np.zeros(shape=(len(t), 1), dtype=complex)
u0[:, 0] = math.sqrt(P0) * 2 ** (-((1 + 1j * C) / 2) * (2 * (t - t0) / FWHM) ** (2 * m))
u0[:, 0] = np.sqrt(abs(u0[:, 0])**2 + P0 / 30 * np.random.randn(u0.shape[0]),
                                   dtype="complex")
nz_step = [10]
list_output = []
list_delta = []
for nz in nz_step:
    # output
    u1 = Channel.ssprop(u0, dt, dz, nz, alpha, betap, gamma)

    title_graph_1 = 'Plot canal: ' + str(nz * dz) + 'Km, alpha = ' \
                    + str(alpha) + '_beta_2_' + str(beta2) + '_gamma_' \
                    + str(gamma) + '.png'

    Plotter.plot_pulse_input_and_output([[t, abs(u0[:, 0]) ** 2, "Input da fibra"],
                                         [t, abs(u1[:, 0]) ** 2, "Output da fibra"]],
                                        graph_title=title_graph_1,
                                        x_graph='s',
                                        y_graph="")
