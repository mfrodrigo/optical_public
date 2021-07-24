% constants - universal and SOA constants
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

% universal constants

global c e k h hbar kT

c = 3e8; % speed light
h = 6.63e-34; % Planck constant (J-s)
k = 1.38e-23; % Boltzmann constant (J/K)
m0 =9.11e-31; % electron rest mass (kg)
e = 1.6e-19; % electron charge (C)
meV = 1e-3*e; 
T = 300; % Temperature (K)

hbar = h/(2*pi); % planck constant/2pi
kT = k*T;

% SOA parameters

global yAs Lc Lt d W confine Kg n1 n2 neq0 eta_in eta_out R1 R2 K0 K1 Arad Brad Anrad Bnrad Caug Dleak...
    me mhh mlh

yAs = 0.892;  % molar fraction of Arsenide in the active region
Lc = 600; % Central active region length (um)
Lt = 100; % Tapered active region length (um)
d = 0.4e-6;  %Active region thickness (m)
W = 0.4e-6; % Central active region width (m)
confine = 0.45; % Optical confinement factor 
Kg = 0.9e-10; % Bandgap shrinkage coefficient (eVm)
n1 = 3.22; % InGaAsP active region refractive index
n2 = 3.167; % InP region refractive index
dn1dn = -1.8e-26; % Differential of active region refractive index with respect to carrier density (m^-3)
neq0 = 3.22; % Equivalent refractive index at zero carrier density (from (29))
dneqdn = n1*confine/(neq0*(2-confine))*dn1dn; % Differential of equivalent refractive index with respect to carrier density (equation (30)) (m^-3)
eta_in = 3; % Input coupling loss (dB)
eta_out = 3; % Output coupling loss (dB)
R1 = 5e-5; % Input facet reflectivity
R2 = 5e-5; % Output facet reflectivity
K0 = 6200; % Carrier independent absorption loss coefficient (m^-1)
K1 = 7500; % Carrier dependent absorption loss coefficient (m^2)
Arad = 1e7; % Linear radiative recombination coefficient (s^-1)
Brad = 5.6e-16; % Bimolecular radiative recombination coefficient (m^3 s^-1)
Anrad = 3.5e8; % Linear nonradiative recombination coefficient due to traps (s^-1)
Bnrad = 0e-16; % Bimolecular nonradiative recombination coefficient (m^3 s^-1)
Caug = 3e-41; % Auger recombination coefficient (m^6s^-1)
Dleak = 0e48; % Leakage recombination coefficient (m^13.5s^-1)

me = m0*0.045; % Effective mass of electron in the CB
mhh = m0*0.46; % Effective mass of a heavy hole in the VB
mlh = m0*0.056; % Effective mass of a light hole in the VB

% SOA parameters derived from above parameters

L = (Lc + Lt)*1e-6; % average length

r1 = sqrt(R1);   % amplitude reflectivities
r2 = sqrt(R2);

eta_in = 10^(-eta_in/10);  % converting from dB to linear quantities
eta_out = 10^(-eta_out/10);

neq = ((n1^2 - n2^2)*confine/(2+confine) + n2^2)^0.5;
delta_Em = h*c/(2*neq*L);  % longitudinal mode energy spacing (assumed carrier independent) (J) 

% parameters used in nilsson function

global mdh nc nv

mdh = (mhh^1.5 + mlh^1.5)^(2/3);
nc = 2*(me*kT/(2*pi*hbar^2))^1.5;
nv = 2*(mdh*kT/(2*pi*hbar^2))^1.5;

% parameters used in gain_coeff function

global k0

k1 = c^2*(h^1.5)/(4*sqrt(2)*(pi^1.5)*n1^2);
k2 = (2*me*mhh/(hbar*(me + mhh)))^1.5;
k0 = k1*k2;

