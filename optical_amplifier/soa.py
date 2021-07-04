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
        t_end = t[-1]
        h = 0.05e-12  # simulation steps
        h0 = self.g0 * self.L
        esat = 10e-12
        pin = np.ones((1, len(t)))
        tf = np.ones((1, len(t)))

        for i in range(len(t)):
            t_aux = t[i]
            H = u[i]
            t = t - t_end * 1e-12 / 2
            tau0 = FWHM / 1.7627
            tauc = 200e-12  # carrier recover time
            tparam = t / tau0
            tparam1 = (t + (h / 2)) / tau0
            tparam2 = (t + h) / tau0
            Pin = sqrt(ein / (2 * tau0)) * (sech(-tparam + 3.52)) + sqrt(ein / (2 * tau0)) * (
                sech(-tparam + 10.57)) + sqrt(ein / (2 * tau0)) * (sech(-tparam + 17.65)) + sqrt(ein / (2 * tau0)) * (
                      sech(-tparam + 45.83))
            Pin2 = sqrt(ein / (2 * tau0)) * (sech(-tparam1 + 3.52)) + sqrt(ein / (2 * tau0)) * (
                sech(-tparam1 + 10.57)) + sqrt(ein / (2 * tau0)) * (sech(-tparam1 + 17.65)) + sqrt(ein / (2 * tau0)) * (
                       sech(-tparam1 + 45.83))
            Pin3 = sqrt(ein / (2 * tau0)) * (sech(-tparam2 + 3.52)) + sqrt(ein / (2 * tau0)) * (
                sech(-tparam2 + 10.57)) + sqrt(ein / (2 * tau0)) * (sech(-tparam2 + 17.65)) + sqrt(ein / (2 * tau0)) * (
                       sech(-tparam2 + 45.83))

            k1 = (h0 - H) / tauc - (abs(Pin) ^ 2) * (exp(H) - 1) / esat
            k2 = (h0 - (H + 0.5 * h * k1)) / tauc - (abs(Pin2) ** 2) * (exp(H + 0.5 * h * k1) - 1) / esat
            k3 = (h0 - (H + 0.5 * h * k2)) / tauc - (abs(Pin2) ** 2) * (exp(H + 0.5 * h * k2) - 1) / esat
            k4 = (h0 - (H + h * k3)) / tauc - (abs(Pin3) ** 2) * (exp(H + h * k3) - 1) / esat
            u[i + 1] = u[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            t[i + 1] = t[i] + h

            pin[i] = Pin
            tf[i] = t

        gain = exp(u)  # gain of SOA
        phi = -0.5 * self.alpha * u  # phase difference os SOA

        Pout = (abs(pin) ** 2) * gain  # output power of SOA

        return Pout