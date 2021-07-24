% SOA_OP - SOA operating point simulation
% Determines output optical spectrum, noise figure and spatial
% distributions using paper 'Wideband Semiconductor Optical Amplifier
% Steady-State Numerical Model' by Michael J. Connelly, IEEE Journal
% Quantum Electron, vol. 37, no. 3, pp.439-447, 2001.
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

clear all
clc
disp('SOA operating point simulation')
constants;  % load in SOA and universal constants
simul_params;  % load in numerical model parameters and vectors used in spectrum calculations

bias = 150e-3;  % bias current in A

% input signal

use_input_signal = 1; % to include the input signal in the model
display_signal_OSA = 1;  % to display the signal on the OSA
if use_input_signal
    Ns = 2; % number of input signals
    Pin_dBm = [-40.0,-20];  % input power in dBm
    Pin = 1e-3*10.^(Pin_dBm/10);  % input signal power in W
    wavelength_s = [1550,1555]*1e-9;  % signal wavelength in m
    Es = h*c./(wavelength_s); % signal energy
    Ain = sqrt(Pin)./sqrt(Es);  % input signal amplitude (square root of input photon rate)
    beta_s = 2*pi*neq*Es/(h*c);  % signal propagation coefficient
end

% Optical spectrum analyser settings

lambda_res = 0.1e-9;  % OSA resolution bandwidth (m)
lambdaOSA0 = 1500e-9;  % OSA start wavelength (m)
lambdaOSA1 = 1600e-9; % OSA end wavelength (m)
power_scale = 'dBm';  % choose power scale ('dBm' or 'linear')
switch power_scale
    case {'dBm'}
        PmindBm = -40; % OSA minimum power (dBm)
        PmaxdBm = -10; % OSA maximum power (dBm)
    case {'linear'}
        Pminlin = 1e-5;  % OSA minimum power (mW)
        Pmaxlin = 0.2e-1; % OSA maximum power (mW)
end

% OSA settings cannot be outside the mode wavelength range

if lambdaOSA0 < lambda0
    lambdaOSA0 = lambda0;
end
if lambdaOSA1 > lambda1
    lambdaOSA1 = lambda1;
end
 
% model initial conditions

weight(1:Nz) = 0.1;  % weighting factor
n(1:Nz) = 1.2e24; % initial guess for carrier density

if use_input_signal
    Asp = zeros(Nz+1,Ns); % forward travelling signal amplitude
    Asm = zeros(Nz+1,Ns); % backward travelling signal amplitude
end

Nsp = zeros(Nz+1,Nm+1);  % forward travelling ASE  (single polarisation)
Nsm = Nsp; % backward travelling ASE (single polarisation)

oldsignQ(1:Nz) = 1;  % assuming a positive sign for initial value of Q (the algorithm will still converge regardless of reasonable initial conditions)

tol = 999;  %

% Start main iteration 

while(tol > tolerence) % calculate coefficients of the travelling-wave equations    
    
    Nsp_old = Nsp;  % storing old values of Nsp and Nsm for use in tolerance calculation
    Nsm_old = Nsm;
    
    % attenuation coefficient - energy independent    
    
    if use_input_signal
        for J = 1:Ns
            alpha_s(:,J) =  K0 + confine*K1*n/1e24;  % signal attenuation coefficient
        end
    end
    
    for J = 1:Nm+1;
        alpha(:,J) = K0 + confine*K1*n/1e24;  % attenuation coefficient
    end
    
% material gain coefficients and spontaneous emission term in equations (37,38)    
    
    for I = 1:Nz
        if use_input_signal
            dummy = gain_coeff(n(I),Es);
            gm_s(I,:) = dummy(1,:);  % material gain coefficient for signal
        end
        dummy = gain_coeff(n(I),E);
        gm(I,:) = dummy(1,:);  % material gain coefficient for ASE
        Rsp(I,:) = dummy(2,:);   % Spontaneous emission term in equations (37,38)
    end

    %    boundary conditions - signal
    
    if use_input_signal    
        Asp(1,:) = (1-r1)*sqrt(eta_in)*Ain + r1*Asm(1,:);
        Asm(Nz+1,:) = r2*Asp(Nz+1,:);
    end
    
    %    boundary conditions - ASE
    
    Nsp(1) = R1*Nsm(1);  % input facet
    Nsm(Nz+1) = R2*Nsp(Nz+1);  % output facet
    
    % solve travelling-wave equations for signal
    
    if use_input_signal
        for I = 2:(Nz+1)
            Asp(I,:) = Asp(I-1,:).*exp((-j*beta_s + 0.5*(confine*gm_s(I-1,:)-alpha_s(I-1,:)))*dz);
        end

        for I = Nz:-1:1
            Asm(I,:) = Asm(I+1,:).*exp((-j*beta_s + 0.5*(confine*gm_s(I,:)-alpha_s(I,:)))*dz);
        end
    end
        
    % solve travelling-wave equations for ASE
   
    for I = 2:(Nz+1)
        Nsp(I,:) = Nsp(I-1,:).*exp((confine*gm(I-1,:)-alpha(I-1))*dz) +...
            + Rsp(I-1,:).*(exp((confine*gm(I-1,:)-alpha(I-1))*dz)-1)./(confine*gm(I-1,:)-alpha(I-1));
    end
    
    for I = Nz:-1:1
        Nsm(I,:) = Nsm(I+1,:).*exp((confine*gm(I,:)-alpha(I))*dz) +...
            Rsp(I,:).*(exp((confine*gm(I,:)-alpha(I))*dz)-1)./(confine*gm(I,:)-alpha(I));
    end
    
    % normalise noise photon rates
    
    Gs = exp(sum((confine*gm-alpha)*dz));    % single pass gain     
    gamma = 4*Gs*sqrt(R1*R2)./(1-sqrt(R1*R2)*Gs).^2;    
    K = 1./sqrt(1+gamma.^2);
    
    % percentage tolerance using ASE, at energies above maximum bandgap only  
    
    Eg = max(egap(n));    
    index = max(find(E <= Eg)) + 1;  
    
    % tolerence calculation uses ASE throughout the SOA
    
    tol = max(50*max(abs(Nsp(2:(Nz+1),index:Nm+1) - Nsp_old(2:(Nz+1),index:Nm+1))./Nsp(2:(Nz+1),index:Nm+1) +...
        abs(Nsm(1:Nz,index:Nm+1) - Nsm_old(1:Nz,index:Nm+1))./Nsm(1:Nz,index:Nm+1)));
    
%    disp(['error = ' num2str(tol) ' %']);  % display current error
    
    % update weights and carrier density 
    
    for I = 1:Nz        
        if use_input_signal
            Q(I) = bias/(e*d*W*L) - ( (Anrad+Arad)*n(I) + (Bnrad+Brad)*n(I)^2 + Caug*n(I)^3 + Dleak*n(I)^5.5) -...
                0.5*confine/(d*W)*sum(gm_s(I,:).*(abs(Asp(I,:)).^2 + abs(Asp(I+1,:)).^2 + abs(Asm(I,:)).^2 + abs(Asm(I+1,:)).^2)) -...
                confine/(d*W)*sum(gm(I,:).*K.*(Nsp(I,:)+Nsp(I+1,:)+Nsm(I,:)+Nsm(I+1,:)));   % RHS of carrier density rate equation
        else
            Q(I) = bias/(e*d*W*L) - ( (Anrad+Arad)*n(I) + (Bnrad+Brad)*n(I)^2 + Caug*n(I)^3 + Dleak*n(I)^5.5) -...
                confine/(d*W)*sum(gm(I,:).*K.*(Nsp(I,:)+Nsp(I+1,:)+Nsm(I,:)+Nsm(I+1,:)));   % RHS of carrier density rate equation
        end
        
        if sign(Q(I)) ~= oldsignQ(I)   % weight adjustment
            weight(I) = weight(I)/2;
        end
        oldsignQ(I) = sign(Q(I));
        
        if Q(I) > 0    % new carrier density
            n(I) = n(I)*(1+weight(I));
        else
            n(I) = n(I)/(1+weight(I));
        end        
    end  
    
end % of tolerence iteration

disp('Numerical algorithm converged')

% Output signal powers

if use_input_signal
    Aout = (1-r2)*sqrt(eta_out)*Asp(Nz+1,:);  % output signal amplitude from output facet
    Pout = Es.*(abs(Aout).^2); % output signal power from output facet
    Aout_back = (1-r1)*sqrt(eta_in)*Asm(1,:); % output signal amplitude from input facet
    Pout_back = Es.*(abs(Aout_back).^2);  % output power from input facet
    disp(['Signal wavelength(s) = ' num2str(wavelength_s/1e-9) ' nm'])
    disp(['Signal input power(s) = ' num2str(Pin_dBm) ' dBm'])
    disp(['Signal gain(s) = ' num2str(10*log10(Pout./Pin)) ' dB'])
end

% total ASE output powers

Nout = 2*eta_out*(1-R2)*sum(K.*Nsp(Nz+1,:).*E);  % output total ASE power
Nout_back = 2*eta_in*(1-R1)*sum(K.*Nsm(1,:).*E);  % output total backward ASE power
disp(['Total output ASE power  = ' num2str(Nout/1e-3) ' mW' '  Total output ASE power (backward)  = ' num2str(Nout_back/1e-3) ' mW'])

% Calculate and display SOA output power spectrum on an OSA

sigmaN = 2*eta_out*(1-R2)*K.*Nsp(Nz+1,:).*E*h/delta_E;  % output ASE power spectral density (Watts/Hz) for energy vector E

sigmaN_res = interp1(E,sigmaN,Eres,'spline'); % output ASE power spectral density (Watts/Hz) at resonant energies Eres
sigmaN_res(find(sigmaN_res <= 0)) = 0;

% calculating material gain coefficient for energy vector with spacing =
% 1/2 longitudinal mode energy difference

clear gm
for I = 1:Nz
    dummy = gain_coeff(n(I),Efine);
    gm(I,:) = dummy(1,:);  % material gain coefficient
end   

% interpolation of attenuation coefficient for energy vector with spacing =
% longitudinal mode energy difference

[X,Y] = meshgrid(z,E);
[X1,Y1] = meshgrid(z,Efine);
alpha = interp2(X,Y,alpha',X1,Y1,'linear',2)';

% single pass gain  (assumed equal for resonant and antiresonant energies

clear Gs
for J = 1:(Km*Nm*2)
    Gs(J) = exp(sum((confine*gm(:,J)-alpha(:,J))*dz));   
end

% resonant and anti-resonant gain combined in single vector

G(1:2:(Km*Nm*2)) = (1-R1)*(1-R2)*Gs(1:2:(Km*Nm*2))./(1-sqrt(R1*R2)*Gs(1:2:(Km*Nm*2))).^2;  % resonant gain  
G(2:2:(Km*Nm*2)) = (1-R1)*(1-R2)*Gs(2:2:(Km*Nm*2))./(1+sqrt(R1*R2)*Gs(2:2:(Km*Nm*2))).^2;  % antiresonant gain 

Gfine2 = interp1(Efine,G,Efine2,'spline','extrap');  % interpolating to finer energy vector

J = 1;
for I = 1:20:length(Efine2)   % integrating amplifier gain profile over a mode spacing
    Gfine2_int(J) = sum(Gfine2(I:I+19))/(20);  % size of Gfine2_int = size(Eres)
    J = J+1;
end

clear sigmaN_fine   % o/p ASE power spectral density (size = size(Efine))
for J = 1:(length(Efine)-1)
    sigmaN_fine(J) = sigmaN_res(floor(J/2)+1)*G(J)/Gfine2_int(floor(J/2)+1);
end

clear sigmaN_fine2  % o/p ASE power spectral density (size = size(Efine2))
for J = 1:(length(Efine2)-1)
    sigmaN_fine2(J) = sigmaN_res(floor(J/20)+1)*Gfine2(J)/Gfine2_int(floor(J/20)+1);
end

disp('SOA operating point simulation complete')
 
% use moving averages to obtain OSA display

wavelength_uniform = lambda0:0.01e-9:(lambda1-0.01e-9);  % uniform wavelength vector with spacing = 0.01 nm
E_uniform = h*c./wavelength_uniform;
G_uniform = interp1(wavelength_fine2,Gfine2,wavelength_uniform); 
G_uniform = G_uniform*eta_in*eta_out;  % fibre-to-fibre gain used in noise figure calculation

sigma_OSA = interp1(wavelength_fine2(1:(length(Efine2)-1)),sigmaN_fine2,wavelength_uniform);  % interpolating power spectral density to uniform wavelength spacing
NF = sigma_OSA ./(E_uniform.*G_uniform) + eta_out./G_uniform;  % noise figure (equation 63)

sigma_OSA = sigma_OSA*c./wavelength_uniform.^2;  % power spectral density (Watts/m)

windowSize = ceil(lambda_res/0.01e-9);  % window for averaging over OSA resolution bandwidth
sigma_OSA = filter(ones(1,windowSize)/windowSize,1,sigma_OSA)*lambda_res;  % power in a wavelength spacing = lambda_res

if use_input_signal & display_signal_OSA
    for J = 1:Ns
        index = min(find(wavelength_uniform >= wavelength_s(J)));
        sigma_OSA(index) = sigma_OSA(index) + Pout(J);  % adding signal to power spectral density for display on OSA
    end
end

sigma_OSA(find(sigma_OSA <= 0)) = 1e-9; % making zero values in the spectrum equal to a small number to avoid problem in taking the log 

% display SOA output spectrum on an OSA

figure(1)
switch power_scale
    case {'dBm'}
        plot(wavelength_uniform/1e-9,db(sigma_OSA/1e-3,'power'))
        xlabel('Wavelength (nm)','Fontsize',14)
        ylabel('Power (dBm)','Fontsize',14)
        title(['Optical spectrum analyser display of SOA output. Resolution bandwidth = ' num2str(lambda_res/1e-9) ' nm'],'Fontsize',14);   
        axis([lambdaOSA0/1e-9 lambdaOSA1/1e-9 PmindBm PmaxdBm ]) 
    case {'linear'}
        plot(wavelength_uniform/1e-9,sigma_OSA/1e-3)
        xlabel('Wavelength (nm)','Fontsize',14)
        ylabel('Power (mW)','Fontsize',14)
        title(['Optical spectrum analyser display of SOA output. Resolution bandwidth = ' num2str(lambda_res/1e-9) ' nm'],'Fontsize',14);   
        axis([lambdaOSA0/1e-9 lambdaOSA1/1e-9 Pminlin Pmaxlin ]) 
end

% plot noise figure spectrum

figure(2)
plot(wavelength_uniform/1e-9,db(NF,'power'))
axis([lambdaOSA0/1e-9 lambdaOSA1/1e-9 0 25])
xlabel('Wavelength (nm)','Fontsize',14);
ylabel('Noise figure (dB)','Fontsize',14)
title('Noise figure spectrum','Fontsize',14)

% display spatial distributions

figure(3)
subplot(2,1,1)
plot(z/1e-6,n/1e24)
xlabel('Distance from input facet (\mum)','Fontsize',14);
ylabel('Carrier density (10^2^4 m^-^3)','Fontsize',14)
title('Carrier density spatial distribution','Fontsize',14)

subplot(2,1,2)
plot(z1/1e-6,sum(Nsp'),':',z1/1e-6,sum(Nsm'),'.',z1/1e-6,sum(abs(Asp)'.^2),z1/1e-6,sum(abs(Asm)'.^2))
legend('Forward travelling ASE', 'Backward travelling ASE','Forward travelling signal', 'Backward travelling signal','Location','North');
xlabel('Distance from input facet (\mum)','Fontsize',14);
ylabel('Total photon rate (s^-^1)','Fontsize',14)
title('Photon rates spatial distributions','Fontsize',14)
axis([0 L/1e-6 0 inf])     