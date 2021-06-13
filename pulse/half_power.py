import numpy as np


def return_half_power(t, u):
    u_squared_module = abs(u[:, 0]) ** 2
    max_peak = np.amax(u_squared_module)
    for i in range(u_squared_module.shape[0]):
        if u_squared_module[i] >= max_peak / 2:
            aux_1 = i
            break
    for j in range(aux_1, u_squared_module.shape[0]):
        if u_squared_module[j] <= max_peak / 2:
            aux_2 = j
            break

    delta_t = t[aux_2] - t[aux_1]
    return [t[aux_1], t[aux_2], delta_t, max_peak, u_squared_module[aux_1], u_squared_module[aux_2]]

