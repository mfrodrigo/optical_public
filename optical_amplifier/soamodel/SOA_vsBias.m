% SOA_vsBias - SOA simulations versus bias current
% Determines signal gain and approximate noise figure versus bias current
% distributions using paper 'Wideband Semiconductor Optical Amplifier
% Steady-State Numerical Model' by Michael J. Connelly, IEEE Journal
% Quantum Electron, vol. 37, no. 3, pp.439-447, 2001.
%   Author: M. Connelly, University of Limerick, Ireland.
%   $Revision: 1.0. $Date: 17/5/2007

clear all
clc
disp('SOA versus bias current simulation')
constants;  % load in SOA and universal constants
simul_params;  % load in numerical model parameters and vectors used in spectrum calculations

bias = [10:20:160]*1e-3;  % bias current vector in A

% input signal

Pin_dBm = -40.0;  % input power in dBm
Pin = 1e-3*10.^(Pin_dBm/10);  % input signal power in W
wavelength_s = 1550*1e-9;  % signal wavelength in m
Es = h*c./(wavelength_s); % signal energy
Ain = sqrt(Pin)./sqrt(Es);  % input signal amplitude (square root of input photon rate)
beta_s = 2*pi*neq*Es/(h*c);  % signal propagation coefficient
 
% model initial conditions

for Ibias = 1:length(bias)
    
    disp(['Bias current = ' num2str(bias(Ibias)/1e-3) ' mA'])
    weight(1:Nz) = 0.1;  % weighting factor
    n(1:Nz) = 1.2e24; % initial guess for carrier density
    
    Asp = zeros(Nz+1); % forward travelling signal amplitude
    Asm = zeros(Nz+1); % backward travelling signal amplitude
    
    Nsp = zeros(Nz+1,Nm+1);  % forward travelling ASE  (single polarisation)
    Nsm = Nsp; % backward travelling ASE (single polarisation)
    
    oldsignQ(1:Nz) = 1;  % assuming a positive sign for initial value of Q (the algorithm will still converge regardless of reasonable initial conditions)
    
    tol = 999;  %
    
    % Start main iteration 
    
    while(tol > tolerence) % calculate coefficients of the travelling-wave equations    
        
        Nsp_old = Nsp;  % storing old values of Nsp and Nsm for use in tolerance calculation
        Nsm_old = Nsm;
        
        % attenuation coefficient - energy independent    
        
        alpha_s =  K0 + confine*K1*n/1e24;  % signal attenuation coefficient
        
        for J = 1:Nm+1;
            alpha(:,J) = K0 + confine*K1*n/1e24;  % attenuation coefficient
        end

        % material gain coefficients and spontaneous emission term in equations (37,38)
        
        for I = 1:Nz
            dummy = gain_coeff(n(I),Es);
            gm_s(I) = dummy(1,:);  % material gain coefficient for signal
            dummy = gain_coeff(n(I),E);
            gm(I,:) = dummy(1,:);  % material gain coefficient for ASE
            Rsp(I,:) = dummy(2,:);   % Spontaneous emission term in equations (37,38)
        end

        %    boundary conditions - signal
        
        Asp(1) = (1-r1)*sqrt(eta_in)*Ain + r1*Asm(1);
        Asm(Nz+1) = r2*Asp(Nz+1);
        
        %    boundary conditions - ASE
        
        Nsp(1) = R1*Nsm(1);  % input facet
        Nsm(Nz+1) = R2*Nsp(Nz+1);  % output facet
        
        % solve travelling-wave equations for signal
        
        for I = 2:(Nz+1)
            Asp(I) = Asp(I-1).*exp((-j*beta_s + 0.5*(confine*gm_s(I-1)-alpha_s(I-1)))*dz);
        end

        for I = Nz:-1:1
            Asm(I) = Asm(I+1).*exp((-j*beta_s + 0.5*(confine*gm_s(I)-alpha_s(I)))*dz);
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
        
%         disp(['error = ' num2str(tol) ' %'])  % display current error
        
        % update weights and carrier density 
        
        for I = 1:Nz        
            Q(I) = bias(Ibias)/(e*d*W*L) - ( (Anrad+Arad)*n(I) + (Bnrad+Brad)*n(I)^2 + Caug*n(I)^3 + Dleak*n(I)^5.5) -...
                0.5*confine/(d*W)*(gm_s(I).*(abs(Asp(I)).^2 + abs(Asp(I+1)).^2 + abs(Asm(I)).^2 + abs(Asm(I+1)).^2)) -...
                confine/(d*W)*sum(gm(I,:).*K.*(Nsp(I,:)+Nsp(I+1,:)+Nsm(I,:)+Nsm(I+1,:)));   % RHS of carrier density rate equation          

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
    
    Aout = (1-r2)*sqrt(eta_out)*Asp(Nz+1);  % output signal amplitude from output facet
    Pout(Ibias) = Es.*(abs(Aout).^2); % output signal power from output facet
   
    % total ASE output powers
    
    Nout(Ibias) = 2*eta_out*(1-R2)*sum(K.*Nsp(Nz+1,:).*E);  % output total ASE power
    
    % Calculate and displace SOA output power spectrum on an OSA
    
    sigmaN_spec = 2*eta_out*(1-R2)*K.*Nsp(Nz+1,:).*E*h/delta_E;  % output ASE power spectral density (Watts/Hz) for energy vector E
    
    sigmaN(Ibias) = interp1(E,sigmaN_spec,Es,'spline'); % interpolating to find output ASE power spectral density at signal wavelength   
    
end

Gain =  Pout/Pin;  % fibre-to-fibre gain
NF = sigmaN./(Es.*Gain) + eta_out./Gain;  % noise figure (equation 63)

disp('SOA versus bias simulation complete')

% plot signal gain and noise figure vs. bias 

figure(1)
[AX,H1,H2] = plotyy(bias/1e-3,db(Gain,'power'),bias/1e-3,db(NF,'power'),'plot');
xlabel('Bias current (mA)','Fontsize',14);
title(['Fibre-to-fibre gain and noise figure at ' num2str(wavelength_s/1e-9) ' nm. Input signal power = ' num2str(Pin_dBm) ' dBm'],'Fontsize',14)

set(get(AX(1),'Ylabel'),'String','Fibre-to-fibre gain (dB)','Fontsize',14)
set(AX(1),'XLim',[0 max(bias)/1e-3]); 
set(AX(1),'YLim',[-10 30],'Ytick',[-10:2:30]); 
set(H1,'LineStyle','-','Marker','+')

set(get(AX(2),'Ylabel'),'String','Noise figure (dB)','Fontsize',14)
set(H2,'LineStyle','-','Marker','o'); 
set(H2,'LineStyle','-.');
set(AX(2),'XLim',[0 max(bias)/1e-3]); 
set(AX(2),'YLim',[10 30],'Ytick',[10:2:30])
