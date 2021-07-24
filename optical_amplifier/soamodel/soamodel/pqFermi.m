%pqFermi  - plot quasi_Fermi levels
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

plot((0.1:0.01:5)*1e24,nilsson((0.1:0.01:5)*1e24,'conduction')/meV,(0.1:0.01:5)*1e24,nilsson((0.1:0.01:5)*1e24,'valence')/meV);
xlabel('Carrier density (m^-^3)');
ylabel('quasi-Fermi level (eV)');
legend('Conduction band', 'Valence band');
title('quasi-Fermi levels');
