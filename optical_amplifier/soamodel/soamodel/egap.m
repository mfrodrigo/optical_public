function y=egap(n);
% EGAP - approximation to band gap energy (J)
% n - carrier n (m^-3)

%   Author: M. Connelly, University of Limerick, Ireland
%   $Revision: 1.0. $Date: 17/5/2007

global e yAs Kg

y = e*(1.35 - 0.775*yAs + 0.149*yAs^2 - Kg*n.^(1/3));