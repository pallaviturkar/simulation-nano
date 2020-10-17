function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
%clear all
% ******************************************
% Load In experimental data
% ******************************************
BF = importdata('BF.xlsx', 1, 0);
BF_RADIUS = 0.00154;
loc = 50;
BF_t = BF.data(1:loc,1);                  % Seconds
BF_theta = BF.data(1:loc,2);              %   Degrees
%dtdt = gradient(BF_theta,t);
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
V0= 0.03;
lams =0.25;
tss = 0.0001;
Rs(1) = BF_RADIUS;
Theta(1) = BF_theta(1)*pi/180;% radians

%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.1;
seconds = floor(BF_t(end)) + 10;
loop = seconds./dt;
interval = 10; %

%%%%%%%%%%%%
% Liquid Parameters
Dg = 1.24E-05;
rhol = 758;
Csat = 0.237;
L =846*1000;
lamff = 0.171;
mwl = 46;    % Molecular weight of the liquid
Tb = 77+273; % Boiling point of the fuel;

% NANO PROPERTIES
mwn = 101.96     % 26.98(kg/kmol)
rhon = 3950  %2700;%(kg/m^3)% bulk
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
numdr = 500;
R = Rs(1);
dr = R./numdr;
r = 0:dr:Rs(1);
count =0;
C = rhon*V0.*ones(1,size(r,2)-1);
u = 0.*C;
RS = R/cos(pi/2-Theta);
h_max = -(sqrt(RS.^2 - R.^2)-RS);
h = RS*cos(asin(r(1:end-1)./RS)) - (RS-h_max);
J = 0.0001*ones(loop,numdr);

% Volume calculation
Beta = (1-cos(Theta)).^2.*(2+cos(Theta));
V(1) = (R/sin(Theta)).^3*pi*(Beta)/3;

% Csat =Volfg.*rhogf(T)
Tls = T0.*ones(1,size(r,2)-1);
% Theta = 39;

for i=  2:loop
    if i == 2
        count=count+1;
        loc(count) = i;
        t(count) = 0;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        Thetasave(count) = Theta;
        Vsave(count) = V;
        V_dotsave(count) = 0;
        usave(count,:)= u;
        Tlssave(count,:)=Tls;
    end
    
    %Evaporation flux
    lambda = 1/2 - Theta/pi;
    J0 = (Dg.*( 1-Yinf).*Csat./R).*(0.27*Theta.^2 + 1.30).*(0.6381 -0.2239.*(Theta - pi/4).^2);
    J = J0.*(1-r(1:end-1).^2./R.^2).^-lambda;
    
    % Volume fraction
    Volf = C./rhon;
    
    % volume flux Vf
    VolFlux = 2*pi*r(1:end-1).*dr.*J./rhol;
    V_dot = -2.*pi.*sum(J.*r(1:end-1).*dr./rhol);
    V =  V + V_dot.*dt;
    if V < 0
        count=count+1;
        loc(count) = i;
        t(count) = i*dt;
        hsave(count,:) = h;
        Thetasave(count) = Theta;
        Vsave(count) = V;
        V_dotsave(count) = V_dot;
        VolFlux_s(count,:) = VolFlux;
        Jsave(count,:) = J;
        % Csatsave(count,:) = Csat;
        Volfsave(count,:)=Volf;
        Tlssave(count,:)=Tls;
        break
    end
     disp("Volume " + V);
    %     fun = @(x)((3*V*(sin(x)/R)^3)/(pi.*((1-cos(x))^2.*(2+cos(x)))));
    %     Theta = fzero(fun,Theta);
    
    c1 = (3*V)/(pi*R^3);
    syms x
    eqn = (((1-cos(x))^2)*(2+cos(x)))/(sin(x))^3 == c1;
    S = solve(eqn, x, 'Real', true);
    Theta = double(S);
    disp("Theta " + Theta);
    %Theta = asin(V./(pi.*R^3))^0.5;
    
    if mod((dt*i),interval) ==0 | i==10
        count=count+1;
        disp(['saved at ' num2str(i*dt) ' seconds'] );
        disp([' mean (Csat) ' num2str(mean(Csat))]);
        loc(count) = i;
        t(count) = i*dt;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        Vsave(count) = V;
        V_dotsave(count) = V_dot;
        VolFlux_s(count,:) = VolFlux;
        Jsave(count,:) = J;
        Thetasave(count) = Theta;
        Volfsave(count,:)=Volf;
    end
end
% Plotting
loc = 1:count;
figure(1)
plot(t(loc),Vsave(loc),'LineWidth',2)
hold on
plot(BF_t,volume,'r')
xlabel('Seconds');
ylabel('Volume (m^3)')
hold on
figure(2)
plot(t(loc),Thetasave(loc)*180/pi,'LineWidth',2)
hold on
plot(BF_t,BF_theta ,'r')
p=polyfit(t(loc),Thetasave(loc)*180/pi,1);
plot(t(loc),polyval(p,t(loc)),'-.k','LineWidth',2);
xlabel('Seconds');
ylabel('Contact Angle (Degrees)')
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
    ylabel(' Radial velocity(g/m^2s)');
    figure(6)
    plot(r(1:end-1)./R,Tlssave(j,:),c)
end
