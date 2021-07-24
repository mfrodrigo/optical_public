function y=nilsson(n,band);
% Nilsson - approximations to quasi-Fermi level relative to band edge
% n - electron or hole n (m^-3)
% band - "conduction" or "valence"
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

global kT mdh nc nv

switch band
    case {'conduction'}
        delta = n/nc;
        y = kT*(log(delta) + delta.*(64 + 0.05524*delta.*(64 + sqrt(delta))).^-0.25);   
    case {'valence'}
        delta = n/nv;
        y = -kT*(log(delta) + delta.*(64 + 0.05524*delta.*(64 + sqrt(delta))).^-0.25);   
    otherwise
        error('Band must be "conduction" or "valence"'.');
end