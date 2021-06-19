# test Channel

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
from Channel.channel import Channel
from pulse.half_power import return_half_power
from output.plotter import Plotter
from output.tables import Tables

# dt
T = 1000  # (ps) deve ser pelo 4x FWHM
num_samplesperbit = 2 ** 11  # should be 2^n
dt = T / num_samplesperbit  # sampling time(ps) # time step (ps)
t = (np.array(range(1, num_samplesperbit + 1)) - (num_samplesperbit + 1) / 2) * dt

# Pulse
FWHM = 100
t0 = 0
C = 0
m = 1
P0 = 1
u0 = np.zeros(shape=(len(t), 1), dtype=complex)
u0[:, 0] = math.sqrt(P0) * 2 ** (-((1 + 1j * C) / 2) * (2 * (t - t0) / FWHM) ** (2 * m))
dz = 0.5  # distance stepsize (km)

wavelength = 1550  # nm
speed_of_light = 299792.458  # nm/ps

# Fiber 1
D = 17  # [ps/nm.km]
beta2 = -(D * wavelength ** 2) / (math.pi * speed_of_light)  # beta2 (ps^2/km)
betap = np.transpose(np.array([0, 0, beta2]).reshape(1, 3))  # dispersion polynomial
gamma = 0.01
alpha = 0.2/4.343

# DCE Fiber
D_DCE = -100  # [ps/nm.km]
beta2_DCE = -(D_DCE * (wavelength ** 2)) / (math.pi * speed_of_light)  # beta2 (ps^2/km)
betap_DCE = np.transpose(np.array([0, 0, beta2_DCE]).reshape(1, 3))  # dispersion polynomial
gamma_DCE = 0.03
alpha_DCE = 0.4/4.343

nz_step = [200]
list_output_1 = []
list_delta_1 = []
list_output_2 = []
list_delta_2 = []

for nz in nz_step:
    # output
    u1 = Channel.ssprop(u0, dt, dz, nz, alpha, betap, gamma)

    list_values = return_half_power(t, u1)
    list_output_1.append(list_values[3])
    list_delta_1.append(list_values[2])

    title_graph_1 = 'Plot canal: ' + str(nz * dz) + 'Km, alpha = ' \
                    + str(alpha) + '_beta_2_' + str(beta2) + '_gamma_' \
                    + str(gamma) + '.png'

    Plotter.plot_pulse_input_and_output(t, u0, u1, title_graph_1)

    nz_DCE = -int(nz * D / D_DCE)
    u2 = Channel.ssprop(u1, dt, dz, nz_DCE, alpha_DCE, betap_DCE, gamma_DCE)

    list_values = return_half_power(t, u2)
    list_output_2.append(list_values[3])
    list_delta_2.append(list_values[2])

    title_graph_2 = 'Plot canal: ' + str(nz_DCE * dz) + 'Km, alpha = ' + str(alpha_DCE) + '_beta_2_' \
                    + str(beta2_DCE) + '_gamma_' + str(gamma_DCE) + '.png'

    Plotter.plot_pulse_input_and_output(t, u1, u2, title_graph_2)

title_graph_1 = 'Plot Valores de Saída: ' + 'beta_2_' + str(beta2) \
                + '_alpha_' + str(alpha) + '_gamma_' + str(gamma) + '.png'
Plotter.plot_power_output_and_delta_output(np.array(nz_step) * dz,
                                           np.array(list_output_1) / P0 * 100,
                                           list_delta_1,
                                           title_graph_1)

title_graph_2 = 'Plot Valores de Saída: ' + 'beta_2_' + str(beta2_DCE) \
                + '_alpha_' + str(alpha_DCE) + '_gamma_' + str(gamma_DCE) + '.png'
Plotter.plot_power_output_and_delta_output(np.array(nz_step) * (-D / D_DCE) * dz,
                                           np.array(list_output_2) / P0 * 100,
                                           list_delta_2,
                                           title_graph_2)

title_table_1 = 'resultados_de_beta_2_' + str(beta2) + '_alpha_' + str(alpha) + '_gamma_' + str(gamma) + '.csv'
Tables.export_table_results(list(np.array(nz_step) * dz),
                            list_output_1,
                            list_delta_1,
                            title_table_1)

title_table_2 = 'resultados_de_beta_2_' + str(beta2_DCE) + '_alpha_' + str(alpha_DCE) + '_gamma_' + str(
    gamma_DCE) + '.csv'
Tables.export_table_results(list(np.array(nz_step) * (-D / D_DCE) * dz),
                            list_output_2,
                            list_delta_2,
                            title_table_2)
