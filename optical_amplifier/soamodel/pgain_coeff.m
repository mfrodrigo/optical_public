% pgain_coeff - plot gain coefficient spectrum
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

density = 1.5e24;

dummy = gain_coeff(density,E);

[AX,H1,H2] = plotyy(wavelength/1e-9,dummy(1,:),wavelength/1e-9,dummy(2,:));
xlabel('Wavelength (nm)','Fontsize',14);
set(get(AX(1),'Ylabel'),'String','g_m - Material gain coefficient (m^-^1)','Fontsize',12)
set(get(AX(2),'Ylabel'),'String','R_s_p - Additive spontaneous emission term  (s^-^1 m^-^1)','Fontsize',12)