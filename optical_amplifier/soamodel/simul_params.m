% SIMUL_PARAMS  - numerical algorithm parameters
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

tolerence = 0.1;  % desired precentage tolerence
Nz = 100; % number of spatial divisions
Nm = 100; % number of spectrum slices
lambda0 = 1300;  % start wavelength for gain coefficient and ASE spectrum (nm)
lambda1 = 1650; % end wavelength for gain coefficient and ASE spectrum (nm)

% derived parameters

dz = L/Nz; % spatial step
z = dz/2:dz:(L-dz/2);  % distance for carrier density
z1 = 0:dz:L; % distance for photon densities

lambda0 = lambda0*1e-9;  % converting to meters
lambda1 = lambda1*1e-9;  % converting to meters

E0 = h*c/lambda1; % lower energy limit for gain coefficient and ASE spectrum
E1 = h*c/lambda0; % upper energy for gain coefficient and ASE spectrum
global delta_E; 
delta_E = (E1-E0)/Nm;  % energy spacing 

Km = ceil(delta_E/delta_Em);  % number of longitudinal modes per model energy spacing

E1 = E0 + Nm*Km*delta_Em; % energy vector with spacing = longitudinal mode energy spacing 

delta_E = (E1-E0)/Nm;  % energy spacing an integer multiple of the longitudinal mode energy spacing

E = E0:delta_E:E1; % energy vector
wavelength = h*c./E; % wavelength vector

% resonant and anti-resonant energy vectors for ASE spectrum calculations
% assuming E0 corresponds to a resonance (longitudinal mode)

Eres = E0:delta_Em:(E1-delta_Em);  
wavelength_res = h*c./Eres;
Eares = Eres + delta_Em/2;  

Efine = E0:delta_Em/2:(E1-delta_Em/2);  % energy vector of resonant and antiresonant energies
wavelength_fine = h*c./Efine;

Efine2 = E0:delta_Em/20:(E1-delta_Em/20);  % energy vector with spacing = 1/10 of longitudinal mode energy spacing
wavelength_fine2 = h*c./Efine2;
