function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
%clear all
%%
% ******************************************
% Load In experimental data
% ******************************************
BF = importdata('BF.xlsx', 1, 0);
%R = BF.data(4:loc,3); %0.0032/2;
BF_RADIUS = 0.00154;
%t = BF.data(4:end,1);% Seconds
loc = 50;
t = BF.data(1:loc,1);                  % Seconds
BF_theta = BF.data(1:loc,2);           %   Degrees
dtdt = gradient(BF_theta,t);
%R = BF.data(4:loc,3)/1000;             % mm to m
height = BF.data(1:loc,4)/1000;         % mm to m
volume = BF.data(1:loc,5)/(1000^3);     %   mm^3 to m^3
surfaceArea = BF.data(4:loc,6)/(1000^2);% mm^2 to m^2
mass = BF.data(1:loc,7);                % kg
dmdt = gradient(mass,t);                % kg/s
%%% Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% System
P =101325;                              % Pa
Yinf =0;                                %Gas Mass Fraction @ infinity
T0 = 273.15 + 25;%; 21.85 65
V0 = 0.03;
tf=68;
% if V0==0
% Rs(1) = RADIUS;                      % - 0.2*BF_RADIUS;
%    Theta(1) = BF_theta(1)*pi/180;       % radians
% else
% Rs(1) = NF_RADIUS;
% Theta(1) = NF_theta(1)*pi/180;          %radians
% end
seconds = floor(t(end)) + 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.001;
loop = seconds./dt;
interval = 10; %
Tinf = 273.15 + 25;
T0=323
lams =0.25;
hs = 0.0015; 
Rs(1) = BF_RADIUS;  
Theta(1) = BF_theta(1)*pi/180;% radians
%1.5 mm thickness
% Gas Parameters
mwg =28.97;                  % kg/kmol
Tbg =77.36;
Cpg =1007;                   %J./kgK
rhog = 1.184;                %kg./m.^3
lamg =0.02551;               %W./mK
visg =1.849e-5 ;             %(Ns./m.^2) viscosity of air
% Liquid Parameters
mwl =46.07;                  % (kg./kmol)
Tbf =351.8;                  % in K
Dg = 1.24E-05;               %(m.^2./s) Diffusion of liquid into air from FSG
Cpf = 1600;                  %(J./kgK) Specific heat of liquid in vapor phase
Cpl =2870;                   %(J./kgK)Specific heat of liquid as liquid
rhol = 758;                  % (kg./m.^3)
rhogf = 0.145;               % from THermal FLuids Online....
Pv = 29330.9;                % (Pa) vapor pressure at T and 1 atm in pascals
L =846*1000;
%(J./kg) Latent Heat fit
Csat = 0.237;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3rd page
lamf = 0.171;      %approx value       %(W./mK) Thermal conductivity of fuel vapor
visf = 8.5e-6;                         %(Ns./m.^2) Viscosity of fuel vapor
visfl = 694e-6;                        % Ns./m.^2 from Thermal Fluids ONline
epsilon =(.58e-9);                     % approx size (max dimension) of liquid molecule in m
lamff = 0.171;                         % liquid conductivity W/mK
s = 0.02197;                           % surface tension N/m
% NANO PROPERTIES
mwn = 26.98;                              % (kg/kmol)
rhon = 2700;                              %(kg/m^3)% bulk
Cpn = (10^3)*0.91;                        %(J/kgK)
volfM = pi()/6;                           % max volume fraction
Rnano =(1/2)*60e-9;                       % radius of avg nano particl e in meters
df = 1.8;
Rg = (1/2)*250e-9;
ra =((Rg/Rnano)^df*Rnano^2)^(1/2);        % from Keblinski Surface Area Equivilent radius
volf =@(T,Y) (rhon*(1/Y-1)/rhol(T)+1)^(-1);
inVolf = (Rg/Rnano)^(df-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4th page
MaxVolf = volfM*inVolf
%rhon = rhon*MaxVolf;% nano
rholMix =@(T,Yn) ((1-Yn)./rhol(T)+Yn./rhon).^-1;
% Simulation
numdr = 500;
R = Rs(1);
dr = R./numdr;
r = 0:dr:Rs(1);
J = 0.0001*ones(loop,numdr);
RS = R/cos(pi/2-Theta);
h_max = -(sqrt(RS.^2 - R.^2)-RS);
%h = RS*cos(asin(r(1:end-1)./RS)) - (RS-h_max);
%Cap = (2*s(T0)/rhol(T0)/9.81)^0.5;
C = rhon*V0.*ones(1,size(r,2)-1); %Concentration from volume fraction
U = 0.*C;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume calculation
Beta = (1-cos(Theta)).^2.*(2+cos(Theta));
V(1) = (R/sin(Theta)).^3*pi*(Beta)/3;
T = T0.*ones(1,size(r,2)-1);
count =0;
for i=2:500
    % SAVE 
    count=count+1;
    if i == 1
        count=count+1; 
        disp(['saved at ' num2str(i*dt) ' seconds'] ); 
        loc(count) = i; 
        t(count) = 0;
        Tlsave(count,:) = T;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        Vsave(count) = V;
        V_dotsave(count) = 0;
%        thetasave(count) = theta;
        Csave(count,:) = C;
        Usave(count,:) = U;
        VNsave(count) = sum(C.*2*pi.*r(1:end-1).*dr.*h);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5th page
         
        Cmodsave(count,:) = 1.*ones(1,size(r,2)-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Evaporation flux
    lambda = 1/2 - Theta/pi;               % found to be best match to FEM results
    J0 = (Dg.*( 1-Yinf).*Csat./R).*(0.27*Theta.^2 + 1.30).*(0.6381 -0.2239.*(Theta - pi/4).^2);
    J = J0.*(1-r(1:end-1).^2./R.^2).^-lambda;
    % volume V
    % volume flux Vf
    VolFlux = 2*pi*r(1:end-1).*dr.*J./rhol;
    V_dot = -2.*pi.*sum(J.*r(1:end-1).*dr./rhol);
    V = V + V_dot.*dt;   
    Vsave(count) = V
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for rloc = 2:size(r,2)-1
    % U(rloc) = (-1/rhol(T())./r(rloc)./h(rloc)).*sum(r(1:rloc).*dr.*(J(1:rloc) +rhol(T(1:rloc)).*delh(1:rloc)./dt));
    % end
    
    % Radial velocity
    tt=i*dt/tf;
    rt=r(1:end-1)./R;
    u=(1/4).*(1./(1-tt)).*1./rt.*((1-rt.^2).^-lambda - (1-rt.^2));
    U1=u.*R./tf;
    U=U1;
    U(1) = 0;
    loc=find(isnan(U));
    if ~isempty(loc)
        %keyboard
        U(loc)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cprev=C;
% for rloc = 2:size(r,2)-1
% C(rloc) = (C(rloc-1)*U(rloc-1)*dt*2*pi*r(rloc-1)*h(rloc-1)- C(rloc) *U(rloc) *dt*2*pi*r(rloc) *h(rloc)+ C(rloc) *dr *2*pi*r(rloc) *h(rloc))/(2*pi*r(rloc)*h(rloc)*dr);
% end
% Cmod = C./Cprev;
% if ~isempty(find(C<0))

% keyboard
% Tn = Tinf - ((L(T,1).*J).*( h./lamff(T) + hs./lams));
 loc=1:count
 plot(t,J(1:length(t)))
 plot(t,Vsave(1:length(t)))
 hold on
 plot(t,volume)
%  
% hold off
%plot(t,u(1:length(t)))
end


