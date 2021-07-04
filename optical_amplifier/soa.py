"""

"""
import numpy as np
from math import sqrt, exp
from mpmath import sech


class SemiconductorOpticalAmplifier:
    """

    """

    def __init__(self, g0, L, alpha, iter=3200):
        """

        Args:
            g0:
            L:
            alpha:
            iter:
        """
        self.g0 = g0
        self.L = L
        self.alpha = alpha

    def amplify_pulse(self, u, t, ein, FWHM):
        tau0 = FWHM / 1.7627
        tauc = 200e-12  # carrier recover time
        t_end = t[-1]
        h = 0.05e-12  # simulation steps
        h0 = self.g0 * self.L
        esat = 10e-12
        Hr = np.zeros(len(t) + 1)
        pin = np.ones(len(t))
        tf = np.ones(len(t))
        tr = np.ones(len(t) + 1)
        Hr[0] = h0
        tr[0] = t[0]
        for i in range(len(t)):
            t_aux = tr[i]
            H = Hr[i]
            t_aux = t_aux - t_end * 1e-12/2
            tparam = t_aux / tau0
            tparam1 = (t_aux + (h / 2)) / tau0
            tparam2 = (t_aux + h) / tau0
            Pin = sqrt(ein / (2 * tau0)) * (sech(-tparam + 3.52)) + sqrt(ein / (2 * tau0)) * (
                sech(-tparam + 10.57)) + sqrt(ein / (2 * tau0)) * (sech(-tparam + 17.65)) + sqrt(ein / (2 * tau0)) * (
                      sech(-tparam + 45.83))
            Pin2 = sqrt(ein / (2 * tau0)) * (sech(-tparam1 + 3.52)) + sqrt(ein / (2 * tau0)) * (
                sech(-tparam1 + 10.57)) + sqrt(ein / (2 * tau0)) * (sech(-tparam1 + 17.65)) + sqrt(ein / (2 * tau0)) * (
                       sech(-tparam1 + 45.83))
            Pin3 = sqrt(ein / (2 * tau0)) * (sech(-tparam2 + 3.52)) + sqrt(ein / (2 * tau0)) * (
                sech(-tparam2 + 10.57)) + sqrt(ein / (2 * tau0)) * (sech(-tparam2 + 17.65)) + sqrt(ein / (2 * tau0)) * (
                       sech(-tparam2 + 45.83))

            k1 = (h0 - H) / tauc - (abs(Pin) ** 2) * (exp(H) - 1) / esat
            k2 = (h0 - (H + 0.5 * h * k1)) / tauc - (abs(Pin2) ** 2) * (exp(H + 0.5 * h * k1) - 1) / esat
            k3 = (h0 - (H + 0.5 * h * k2)) / tauc - (abs(Pin2) ** 2) * (exp(H + 0.5 * h * k2) - 1) / esat
            k4 = (h0 - (H + h * k3)) / tauc - (abs(Pin3) ** 2) * (exp(H + h * k3) - 1) / esat
            Hr[i + 1] = Hr[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            tr[i + 1] = tr[i] + h

            pin[i] = Pin
            tf[i] = t_aux

        gain = np.exp(Hr)  # gain of SOA
        phi = -0.5 * self.alpha * Hr  # phase difference os SOA

        Pout = (abs(pin) ** 2) * gain[0:len(t)]  # output power of SOA

        return Pout, pin, gain, phi, tf
