function y=gain_coeff(n,E);
% GAIN_COEFF - material and 'spontaneous' gain coefficients
% n - electron or hole density (m^-3)
% E - energy (J)
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

global h kT Arad Brad mhh me k0 confine delta_E

lifetime = (Arad + Brad*n).^-1;  % approximate carrier radiative lifetime

Eg = egap(n);
Ea = (E-egap(n))*mhh/(me + mhh);
Eb = -(E-egap(n))*me/(me + mhh);

Efc = nilsson(n,'conduction');  % quasi-Fermi levels
Efv = nilsson(n,'valence');

fc = (1 + exp((Ea-Efc)/kT)).^-1;  % Fermi-Dirac functions
fv = (1 + exp((Eb-Efv)/kT)).^-1;

for I = 1:length(E)
    if E(I)-Eg > 0
        gm1(I) = sqrt(E(I) - Eg).*fc(I).*(1-fv(I))/(E(I).^2);  % 'spontaneous' gain coefficient 
        gm11(I) = sqrt(E(I) - Eg).*(1-fc(I)).*fv(I)/(E(I).^2);  % absorption coefficient
    else  % for energy less than bandgap
        gm1(I) = 0;
        gm11(I) = 0;
    end
end
gm1 = k0*gm1/lifetime;
gm11 = k0*gm11/lifetime;

gm = gm1 - gm11;  % material gain coefficient
Rsp = confine*gm1*delta_E/h;  % additive spontaneous emission term

y = [gm;Rsp];

