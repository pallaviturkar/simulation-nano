function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
%clear all
% ******************************************
% Load In experimental data
% ******************************************
NanoFluid = importdata('NanoFluid.xlsx', 1, 0);
NF_RADIUS = 0.0024/2;
loc = 42;
NF_t = NanoFluid.data(1:loc,1);                  % Seconds
NF_theta = NanoFluid.data(1:loc,2);
%dtdt = gradient(BF_theta,t);
% height = BF.data(1:loc,4)/1000;         % mm to m
% volume = BF.data(1:loc,5)/(1000^3);     %   mm^3 to m^3
% surfaceArea = BF.data(4:loc,6)/(1000^2);% mm^2 to m^2

%%%%%%%%%%%%%%%%%
% System
P =101325;                              % Pa
Yinf =0;                                % Gas Mass Fraction @ infinity
T0 = 273.15 + 25;%; 21.85 65
Ru =8.314;
tf=40;
Tinf = 273.15 + 25;
V0= 0.02;
lams =0.25;
tss = 0.0001;
Rs(1) = NF_RADIUS;
Theta(1) = NF_theta(1)*pi/180;% radians

%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;
seconds = floor(NF_t(end));
loop = seconds./dt;
interval = 10; %

%%%%%%%%%%%%
% Liquid Parameters
Dg = 1.24E-05;
rhol = 789;
%Csat = 0.39 %%0.237;
L =846*1000;
lamff = 0.171;
mwl = 0.046;    % Molecular weight of the liquid
Tb = 78+273; % Boiling point of the fuel;

% NANO PROPERTIES
mwn = 101.96     % 26.98(kg/kmol)
rhon = 2700  %2700;%(kg/m^3)% bulk
Cpn = 451 %   (10^3)*0.91;%(J/kgK)
volfM = pi()/6;% max volume fraction
Rnano =(1/2)*60e-9;% radius of avg nano particle in meters
df = 1.8;
Rg = (1/2)*250e-9;
% ra =((Rg/Rnano)^df*Rnano^2)^(1/2);% from Keblinski Surface Area Equivilent radius
inVolf = (Rg/Rnano)^(df-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4th page
MaxVolf = volfM*inVolf
%rhon = rhon*MaxVolf;% nano
%rholMix =@(T,Yn) ((1-Yn)./rhol(T)+Yn./rhon).^-1;

% Simulation
numdr = 550;
R = Rs(1);
dr = R./numdr;
r = 0:dr:Rs(1);
C = rhon*V0.*ones(1,size(r,2)-1);
u = 0.*C;
RS = R/cos(pi/2-Theta);
h_max = -(sqrt(RS.^2 - R.^2)-RS);
h = RS*cos(asin(r(1:end-1)./RS)) - (RS-h_max);
%J = 0.0001*ones(loop,numdr);
count =0;

% Volume calculation
Beta = (1-cos(Theta)).^2.*(2+cos(Theta));
V(1) = (R/sin(Theta)).^3*pi*(Beta)/3;

% Csat =Volfg.*rhogf(T)
Tls = T0.*ones(1,size(r,2)-1);
% Ts = T0.*ones(1,size(r,2)-1);
volf =0.02;

for j=  2:loop
    if j == 2
        count=count+1;
        loc(count) = j;
        t(count) = 0;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        Thetasave(count) = Theta;
        Vsave(count) = V;
        V_dotsave(count) = 0;
        %       usave(count,:)= u;
        Tlssave(count,:)=Tls;
        Csave(count,:) = C;
        %       Tssave(count,:)=Ts;
    end
    
    
    %  volf=0;
    % Vapor concentration
    Csat =P./(Ru.*Tls./mwl).*(1-volf).*exp((9.2*100000./(Ru./mwl)).*((1./Tb)-(1./Tls)));
    
    %Evaporation flux
    lambda = 1/2 - Theta/pi;
    J0 = (Dg.*( 1-Yinf).*Csat./R).*(0.27*Theta.^2 + 1.30).*(0.6381 -0.2239.*(Theta - pi/4).^2);
    J = J0.*(1-r(1:end-1).^2./R.^2).^-lambda;
    
    % volume flux Vf
    VolFlux = 2*pi*r(1:end-1).*dr.*J./rhol;
    V_dot = -2.*pi.*sum(J.*r(1:end-1).*dr./rhol);
    V =  V + V_dot.*dt;
    if V < 0
        count=count+1;
        loc(count) = j;
        t(count) = j*dt;
        hsave(count,:) = h;
        Thetasave(count) = Theta;
        Vsave(count) = V;
        V_dotsave(count) = V_dot;
        VolFlux_s(count,:) = VolFlux;
        Jsave(count,:) = J;
        Csatsave(count,:) = Csat;
        volfsave(count,:)=volf;
        Tlssave(count,:)=Tls;
        Csave(count,:) = C;
        break
    end
    
    % Calculation of theta and h at new volume
    c1 = (3*V)/(pi*R^3);
    syms x
    eqn = (((1-cos(x))^2)*(2+cos(x)))/(sin(x))^3 == c1;
    S = solve(eqn, x, 'Real', true);
    Theta = double(S);
    disp("Theta " + Theta);
    
    RS = R/cos((pi/2)-Theta);
    disp("RS " +RS);
    int_h = -(sqrt(RS^2 - R^2)-RS);
    hprev = h;
    h = RS*cos(asin(r(1:end-1)./RS)) - (RS-int_h);
    delh = h - hprev;
    
    % Temperature
    Tls = Tinf - ((L.*J).*( h./lamff + tss./lams));
    Ts  = Tinf- ((L.*J).*( (h+tss)./lams));
    
    for rloc = 2:size(r,2)
        U(rloc) = (-1/rhol./r(rloc)./h(rloc)).*sum(r(2:rloc).*dr.*(J(2:rloc) +rhol.*delh(2:rloc)./dt));
    end
    
    % Radial velocity
    tt=j*dt/tf;
    rt=r(1:end-1)./R;
    u=(1/4).*(1./(1-tt)).*1./rt.*((1-rt.^2).^-lambda - (1-rt.^2));
%     U1=u.*R./tf;
%     U=U1;
    U(1) = 0;
    
    % Particle Concentration
    hsave(j,:) = h;
    for rloc = 2:size(r,2)-1
        Csave(j,rloc) = (dt*2*pi*Csave(j-1,rloc-1)*U(rloc-1)*r(rloc-1)*hsave(j-1,rloc-1)- dt*2*pi*Csave(j-1,rloc) *U(rloc) *r(rloc) *hsave(j-1,rloc)+ dr *2*pi*Csave(j-1,rloc)*r(rloc) *hsave(j-1,rloc))/(2*pi*r(rloc)*hsave(j,rloc)*dr);
        if rloc == 2
            Csave(j, 1) = Csave(j, 2);
        end
    end
    
    %  Volume fraction
    volf = Csave(j,:)./rhon;
   if mod((dt*j),interval) ==0 | i==10
    count=count+1;
    disp(['saved at ' num2str(j*dt) ' seconds'] );
    disp([' mean (Csat) ' num2str(mean(Csat))]);
    loc(count) = j;
    t(j) = j*dt;
    Rsave(j,:) = R;
    Vsave(j) = V;
    V_dotsave(j) = V_dot;
    VolFlux_s(j,:) = VolFlux;
    Jsave(j,:) = J;
    Thetasave(j) = Theta;
    Csatsave(j,:) = Csat;
    usave(j,:)= u;
    volfsave(count,:)=volf;
    Tlssave(j,:)=Tls;
    Tssave(count,:)=Ts;
   end
end
% Plotting
loc = 1:count;

for j =1:count
    c='';
    figure(3)
    plot(r(1:end-1)./R,hsave(j,:)*1000,c)
    xlabel('Normalized Radius');
    ylabel('Height (m)');
    hold on
    figure(4)
    plot(r(1:end-1)./R,Jsave(j,:)*1000,c)
    xlabel('Normalized Radius');
    ylabel('Evaporative Flux J (g/m^2s)');
    hold on
    figure(5)
    plot(r(1:end-1)./R,usave(j,:),c);
    xlabel('Normalized Radius');
    ylabel(' Radial velocity');
    hold on
    figure(6)
    plot(r(1:end-1)./R,Tlssave(j,:),c);
    hold on
    figure(7)
    plot(r(1:end-1)./R,Tssave(j,:),c)
    hold on
    figure(8)
    plot(r(1:end-1)./R,Csave(j,:),c)
    hold on
    figure(9)
    plot(r(1:end-1)./R,volfsave(j,:),c)
    hold on
end
