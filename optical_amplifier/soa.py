"""

"""
import numpy as np
from math import sqrt, exp, pi, log
from mpmath import sech


class SemiconductorOpticalAmplifier:
    """

    """

    c = 3e8  # speed light
    h = 6.63e-34  # Planck constant (J-s)
    k = 1.38e-23  # Boltzmann constant (J/K)
    m0 = 9.11e-31  # electron rest mass (kg)
    e = 1.6e-19  # electron charge (C)
    meV = 1e-3 * e
    T = 300  # Temperature (K)

    hbar = h / (2 * pi)  # planck constant/2pi
    kT = k * T

    yAs = 0.892  # molar fraction of Arsenide in the active region
    Lc = 600  # Central active region length (um)
    Lt = 100  # Tapered active region length (um)
    d = 0.4e-6  # Active region thickness (m)
    W = 0.4e-6  # Central active region width (m)
    confine = 0.45  # Optical confinement factor
    Kg = 0.9e-10  # Bandgap shrinkage coefficient (eVm)
    n1 = 3.22  # InGaAsP active region refractive index
    n2 = 3.167  # InP region refractive index
    dn1dn = -1.8e-26  # Differential of active region refractive index with respect to carrier density (m^-3)
    neq0 = 3.22  # Equivalent refractive index at zero carrier density (from (29))
    dneqdn = n1 * confine / (neq0 * (
            2 - confine)) * dn1dn  # Differential of equivalent refractive index with respect to carrier density (equation (30)) (m^-3)
    eta_in = 3  # Input coupling loss (dB)
    eta_out = 3  # Output coupling loss (dB)
    R1 = 5e-5  # Input facet reflectivity
    R2 = 5e-5  # Output facet reflectivity
    K0 = 6200  # Carrier independent absorption loss coefficient (m^-1)
    K1 = 7500  # Carrier dependent absorption loss coefficient (m^2)
    Arad = 1e7  # Linear radiative recombination coefficient (s^-1)
    Brad = 5.6e-16  # Bimolecular radiative recombination coefficient (m^3 s^-1)
    Anrad = 3.5e8  # Linear nonradiative recombination coefficient due to traps (s^-1)
    Bnrad = 0e-16  # Bimolecular nonradiative recombination coefficient (m^3 s^-1)
    Caug = 3e-41  # Auger recombination coefficient (m^6s^-1)
    Dleak = 0e48  # Leakage recombination coefficient (m^13.5s^-1)

    me = m0 * 0.045  # Effective mass of electron in the CB
    mhh = m0 * 0.46  # Effective mass of a heavy hole in the VB
    mlh = m0 * 0.056  # Effective mass of a light hole in the VB

    # parameters used in nilsson function

    mdh = (mhh ** 1.5 + mlh ** 1.5) ** (2 / 3)
    nc = 2 * (me * kT / (2 * pi * hbar ** 2)) ** 1.5
    nv = 2 * (mdh * kT / (2 * pi * hbar ** 2)) ** 1.5

    def __init__(self):
        pass

    def energy_gap(self, carrier_density):
        """
        this function calculate the approximation to band gap energy (J).

        Returns:

        """
        return self.e * (1.35 - 0.775 * self.yAs + 0.149 * self.yAs ** 2 - self.Kg * carrier_density**(1 / 3))

    def approximation_to_quasi_fermi_level(self, carrier_density, band):
        """

        Args:
            carrier_density:
            band:

        Returns:

        """
        if band == "conduction":
            delta = carrier_density / self.nc
        elif band == "valence":
            delta = carrier_density / self.nv
        else:
            raise NameError('Band must be "conduction" or "valence"')

        return self.kT * (np.log(delta) + delta * (64 + 0.05524 * delta * (64 + np.sqrt(delta))) ** -0.25)

    def gain_coefficient(self, carrier_density, Energy):
        """

        Returns:

        """
        lafetime = (self.Arad + self.Brad * carrier_density) ** -1
        energy_gap = self.energy_gap(carrier_density)
        energy_a = (Energy - energy_gap(carrier_density)) * \
                   self.mhh / (self.me + self.mhh)
        energy_b = (Energy - energy_gap(carrier_density)) * \
                   self.me / (self.me + self.mhh)

        energy_fermi_conduction = \
            self.approximation_to_quasi_fermi_level(carrier_density=carrier_density,
                                                    band='conduction')

        energy_fermi_valence = self.approximation_to_quasi_fermi_level(carrier_density=carrier_density,
                                                                       band='valence')

        fermi_dirac_conduction = (1 + exp(energy_a - energy_fermi_conduction) /
                                  self.kT) ** -1
        fermi_dirac_valence = (1 + exp(energy_b - energy_fermi_valence) /
                               self.kT) ** -1

        gain_coefficient = []
        absorption_coefficient = []
        # for i in range(len(Energy)):
        #     if Energy(i) - energy_gap:
        #         gain_coefficient.append(
        #             sqrt(Energy(i) - energy_gap) *\
        #             fermi_dirac_conduction(i) * (1 - fermi_dirac_valence(i))\
        #             / (Energy(i) ** 2))
        #         absorption_coefficient.append(
        #             sqrt(Energy(i) - energy_gap)*\
        #             (1-fermi_dirac_conduction)*(fermi_dirac_valence(i))
        #
        #         )