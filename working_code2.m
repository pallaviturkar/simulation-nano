function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
%clear all
% ******************************************
% Load In experimental data
% ******************************************
BF = importdata('BF.xlsx', 1, 0);
BF_RADIUS = 0.00154;
loc = 50;
t = BF.data(1:loc,1);                  % Seconds
BF_theta = BF.data(1:loc,2);              %   Degrees
dtdt = gradient(BF_theta,t);
height = BF.data(1:loc,4)/1000;         % mm to m
volume = BF.data(1:loc,5)/(1000^3);     %   mm^3 to m^3
surfaceArea = BF.data(4:loc,6)/(1000^2);% mm^2 to m^2
%%%%%%%%%%%%%%%%%
% System
P =101325;                              % Pa
Yinf =0;                                %Gas Mass Fraction @ infinity
T0 = 273.15 + 50;%; 21.85 65
tf=68;
Tinf = 273.15 + 25;
V0= 0;
lams =0.25;
tss = 0.0001;
Rs(1) = BF_RADIUS;
Theta(1) = BF_theta(1)*pi/180;% radians
%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.1;
seconds = floor(t(end)) + 10;
loop = seconds./dt;
interval = 10; %
%%%%%%%%%%%%
% Liquid Parameters
Dg = 1.24E-05;
rhol = 758;
Csat = 0.237;
L =846*1000;
lamff = 0.171;                         % liquid conductivity W/mK
% NANO PROPERTIES
rhon = 2700;
volfM = pi()/6;
% Simulation
numdr = 500;
R = Rs(1);
dr = R./numdr;
r = 0:dr:Rs(1);
count =0;
C = rhon*V0.*ones(1,size(r,2)-1);
RS = R/cos(pi/2-Theta);
h_max = -(sqrt(RS.^2 - R.^2)-RS);
h = RS*cos(asin(r(1:end-1)./RS)) - (RS-h_max);
% Volume calculation
Beta = (1-cos(Theta)).^2.*(2+cos(Theta));
V(1) = (R/sin(Theta)).^3*pi*(Beta)/3;

for i=  2:500
    if i == 2
        count=count+1;
        loc(count) = i;
        t(count) = 0;
        %        Tlsave(count,:) = T; 
        
        Rsave(count,:) = R;
        hsave(count,:) = h;
        Thetasave(count) = Theta;
        Vsave(count) = V;
        V_dotsave(count) = 0;
    end
    %        for rloc = 2:size(r,2)-1
    %        count=count+1;
    h = RS*cos(asin(r(1:end-1)./RS)) - (RS-h_max);
    %Evaporation flux
    lambda = 1/2 - Theta/pi;               % found to be best match to FEM results
    J0 = (Dg.*( 1-Yinf).*Csat./R).*(0.27*Theta.^2 + 1.30).*(0.6381 -0.2239.*(Theta - pi/4).^2);
    J = J0.*(1-r(1:end-1).^2./R.^2).^-lambda;
    % volume flux Vf
    VolFlux = 2*pi*r(1:end-1).*dr.*J./rhol;
    V_dot = -2.*pi.*sum(J.*r(1:end-1).*dr./rhol);
    V =  V + V_dot.*dt;
    Vsave(count) = V
    if V < 0
        %             disp(['Complete evap at ' num2str(i*dt) 'seconds of simulation time.'])
        count=count+1;
        %             disp(['saved at ' num2str(i*dt) ' seconds'] );
        loc(count) = i;
        t(count) = i*dt;
        % %            Tlsave(count,:) = T;
        %             Rsave(count,:) = R;
        hsave(count,:) = h;
        Thetasave(count) = Theta;
        %  %           hh(count) = int_h;
        %             Vsave(count) = V;
        %             V_dotsave(count) = V_dot;
        %             VolFlux_s(count,:) = VolFlux;
        Jsave(count,:) = J;
        %             Csatsave(count,:) = Csat;
        break
    end
    Theta = asin(V./(pi.*R^3))^0.5;
    
    % Temperature
    Tn = Tinf - ((L.*J).*( h./lamff + tss./lams));
    % Radial velocity
    tt=i*dt/tf;
    rt=r(1:end-1)./R;
    u=(1/4).*(1./(1-tt)).*1./rt.*((1-rt.^2).^-lambda - (1-rt.^2));
    % C(rloc) = (C(rloc-1)*u(rloc-1)*dt*2*pi*r(rloc-1)*h(rloc-1)- C(rloc) *u(rloc) *dt*2*pi*r(rloc) *h(rloc)+ C(rloc));
    % SAVE VALUES
    % ******************************************
    if mod((dt*i),interval) ==0 | i==10
        count=count+1;
        disp(['saved at ' num2str(i*dt) ' seconds'] );
        disp([' mean (Csat) ' num2str(mean(Csat))]);
        %         disp([' mean (tCsat) ' num2str(mean(tCsat))]);
        %         disp([' mean (t2Csat) ' num2str(mean(t2Csat))]);
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         disp([' mean (t3Csat) ' num2str(mean(t3Csat))]);
        %         disp([' mean (Ygf) ' num2str(mean(Ygf))]);
        %         disp([' mean (tYgf) ' num2str(mean(tYgf))]);
        %         disp([' mean (Xgf) ' num2str(mean(Xgf))]);
        %         disp([' mean (tXgf) ' num2str(mean(tXgf))]);
        %         disp([' mean (Xlf) ' num2str(mean(Xlf))]);
        %         disp([' Partial pressure = Xgf (using Pv)' num2str(mean(Pv(T)./P))]);
        loc(count) = i;
        t(count) = i*dt;
        % Tlsave(count,:) = T;
        % Tinf_forcesave(count,:) =Tinf_forced;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        %        hh(count) = int_h;
        Vsave(count) = V;
        V_dotsave(count) = V_dot;
        VolFlux_s(count,:) = VolFlux;
        Jsave(count,:) = J;
        Thetasave(count) = Theta;
        Csatsave(count,:) = Csat;
        % Ygfsave(count,:) = Ygf;
        %         Thetasave(count) = Theta;
        %         MaranSave(count) = Maran;
        %         GravSave(count,:) = Grav;
        %         Csave(count,:) = C;
        %         Xlfsave(count,:) = Xlf;
        %         Ysave(count,:) = Y;
        %         Usave(count,:) = U;
        %         U1save(count,:) = U1;
        %         VNsave(count) = sum(C.*2*pi.*r(1:end-1).*dr.*h);
        %         Cmodsave(count,:) = Cmod;
        %         BL_s(count,:,:) = BL;
        %     end
    end
end
