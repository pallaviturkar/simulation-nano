function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
clear all
%%
% ******************************************
% Load In Data from Sefiane and Bennacer
% @ 25 Celcius
% ******************************************
BaseFluid = importdata('Base_Fluid.xlsx', 1, 0);
BF_RADIUS = 0.0024/2;
BF_t = BaseFluid.data(4:end,1);% Seconds
loc = 50;
BF_t = BaseFluid.data(4:loc,1);% Seconds
BF_theta = BaseFluid.data(4:loc,2);%Degrees
BF_dtdt = gradient(BF_theta,BF_t);
BF_Dradius = BaseFluid.data(4:loc,3)/1000;% mm  to m
BF_height = BaseFluid.data(4:loc,4)/1000;% mm to m
BF_volume = BaseFluid.data(4:loc,5)/(1000^3);%   mm^3 to m^3
BF_surfaceArea = BaseFluid.data(4:loc,6)/(1000^2);% mm^2 to m'^2
BF_mass = BaseFluid.data(4:loc,7);% kg
BF_dmdt = gradient(BF_mass,BF_t);% kg/s
NanoFluid = importdata('Nano_Fluid.xlsx');
NF_RADIUS = 0.0024/2;
NF_t = NanoFluid.data(4:end,1);% Seconds
loc = 50;
NF_t = NanoFluid.data(4:loc,1);% Seconds
NF_theta = NanoFluid.data(4:loc,2);%Degrees
NF_dtdt = gradient(NF_theta,NF_t);
NF_Dradius = NanoFluid.data(4:loc,3)/1000;% mm to m
NF_height = NanoFluid.data(4:loc,4)/1000;% mm to m
NF_volume = NanoFluid.data(4:loc,5)/(1000^3);% mm^3 to m^3
NF_surfaceArea = NanoFluid.data(4:loc,6)/(1000^2);% mm^2 to m^2
NF_mass = NanoFluid.data(4:loc,7);% kg
NF_dmdt = gradient(NF_mass,BF_t);% kg/s
% Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% System
P =101325; % Pa
Yinf =0;%Gas Mass Fraction @ infinity
T0 = 273.15 + 50;%; 21.85 65
V0 = 0.03;
if V0==0
    Rs(1) = BF_RADIUS;% - 0.2*BF_RADIUS;
    Theta(1) = BF_theta(1)*pi/180;% radians
else
    Rs(1) = NF_RADIUS;
    Theta(1) = NF_theta(1)*pi/180;% radians
end
seconds = floor(BF_t(end)) + 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd page
% %%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% System
P =101325; % Pa
Yinf =0;%Gas Mass Fraction @ infinity
T0 = 273.15 + 50;%; 21.85 65
V0 = 0.03;
if V0==0
    
    Rs(1) = BF_RADIUS;% - 0.2*BF_RADIUS;
    Theta(1) = BF_theta(1)*pi/180;% radians
else
    Rs(1) = NF_RADIUS;
    Theta(1) = NF_theta(1)*pi/180;% radians
end
seconds = floor(BF_t(end)) + 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd page

dt = 0.001;
loop = seconds./dt;
interval = 10; %
Tinf = 273.15 + 25;
T0=323;  %%%%%My input
% Solid Surface parameters
% if strcmpi(stringSublayer,'Al')
% lams=237;
% elseif strcmpi(stringSublayer,'Steel')
% lams=50;
% elseif strcmpi(stringSublayer,'PTFE_Mod')
lams = 0.25;%.91;%0.25;%0.25;%237 ;% 237; % conductivity of
% aluminum from in W./mK
% elseif strcmpi(stringSublayer,'PTFE');
% lams =0.25;
% end
%0.91 matches well to base fluid results
hs = 0.0015;%0.0015;%1.5 mm thickness
% Gas Parameters
mwg =28.97;% kg/kmol
Tbg =77.36;
Cpg =@(T) (10.^3).*1.005; %J./kgK
X =[150;200;250;300;350];
Y =[2.3364;1.7458;1.3947;1.1614;.995];
rhog =fit(X,Y,'cubicspline');%kg./m.^3
lamg =@(T) (7.071e-5).*(T-273.15)+2.428e-2;%W./mK
visg =@(T) (4.7225e-8).*(T-273.15)+(1.7238e-5); %(Ns./m.^2) viscosity of air
% Liquid Parameters
mwl =46.07;% (kg./kmol)
Tbf =351.8;% in K
Dg =@(T) (1./100).*(1./100).*(((10.^-3).*T.^1.75.*(1./mwg+1./mwl).^(1./2))./((P./101325).*(50.36.^(1./3)+20.1.^(1./3)).^2)); %(m.^2./s) Diffusion of liquid into air from FSG
Cpf =@(T) (1e3./mwl).*((-8.28925e-5).*T.^2+(0.216104).*T+8.28126);%(J./kgK) Specific heat of liquid in vapor phase
Cpl =@(T) (1e3./mwl).*(98.39+0.5368.*(T-273.15));%(J./kgK)Specific heat of liquid as liquid
rhol =@(T) (-0.8544).*(T-273.15)+806.43;% (kg./m.^3)
rhogf =@(T) exp((-3.3681)+(5.2492e-2).*(T-273.15)+(5.1630e-5).*(T-273.15).^2+(-1.9542e-6).*(T-273.15).^3+(8.6893e-9).*(T-273.15).^4+(-1.1451e-11).*(T-273.15).^5); % from THermal FLuids Online....
Pv =@(T) (133.3224).*((4.0325e-4).*(T-273.15).^3+(2.7952e-2).*(T-273.15).^2+(0.81796).*(T-273.15)+(11.574));% (Pa) vapor pressure at T and 1 atm in pascals
%L=@(T,p)max(1,(.0839.*log(Yp)+1.7831)).*(1e3).*(1e3./mwl).*(50.43).*exp(0.4475.*T./513.9).*(1-T./513.9).^(0.4989);%(J./kg) Latent Heat fit
L =@(T,Xlf)Xlf.*(1e3).*(1e3./mwl).*(50.43).*exp(0.4475.*T./513.9).*(1-T./513.9).^(0.4989);%(J./kg) Latent Heat fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3rd page

lamf =@(T) (4.1841e-7).*T.^2+(-1.6423e-4).*T+.026248;%(W./mK) Thermal conductivity of fuel vapor
visf =@(T) (1.4991e-7)+(3.0741e-8).*T+(-4.4479e-12).*T.^2;
%(Ns./m.^2) Viscosity of fuel vapor
visfl =@(T) (1./1000).*exp((5.8942e-1)+(-2.254e-2).*(T-273.15)+(1.0283e-4).*(T-273.15).^2+(-8.8574e-7).*(T-273.15).^3+(4.7884e-9).*(T-273.15).^4+(-9.7493e-12).*(T-273.15).^5);%Ns./m.^2 from Thermal Fluids ONline
epsilon =(.58e-9); % approx size (max dimension) of liquid molecule in m
lamff =@(T) exp(-1.6976 + -1.2505e-3.*(T-273.15) + 7.59291e-7.*(T-273.15).^2 +5.2361e-8.*(T-273.15).^3 + -3.4986e-10.*(T-273.15).^-4 + 6.4599e-13.*(T-273.15).^5);%liquid conductivity W/mK
s =@(T) 24.419 + -8.1477e-2.*(T-273.15) + -1.1450e-4.*(T-273.1).^2 + 8.6540e-7.*(T-273.15).^3 + -7.6432e-9.*(T-273.15).^4 +1.9148e-11.*(T-273.15).^5;% surface tension N/m
% Mixture Properties of Gas Phase
Sg =1.5.*Tbg;
Sf =1.5.*Tbf;
Sfg =.73.*(Sg.*Sf).^.5;
Sgf =Sfg;
Sff =.73.*(Sf.*Sf).^.5;
Sgg =.73.*(Sg.*Sg).^.5;
A =@(T,vi,vj,Mi,Mj,Si,Sj,Sij).25.*(1+((((vi)./(vj)).*((Mj)./(Mi)).^.75))); %....*((T+Si)./(T+Sj))).^(1./2)).^2).*(T+Sij)./(T+Si);
CpMix =@(T,Ygf) Ygf.*Cpf(T)+(1-Ygf).*Cpg(T); % linear mixture rule
rhogMix =@(T,Ygf) Ygf.*rhogf(T)+(1-Ygf).*rhog(T)'; % linear mixture rule
% Gas Lambda Mix Calc Variables ------------------------------------------
Sg=1.5*Tbg;
Sf=1.5*Tbf;
Sfg=.73*(Sg*Sf)^.5;
Sgf=Sfg;
Sff=.73*(Sf*Sf)^.5;
Sgg=.73*(Sg*Sg)^.5;
A=@(T,vi,vj,Mi,Mj,Si,Sj,Sij).25*(1+((((vi)/(vj))*((Mj)/(Mi))^.75*((T+Si)/(T+Sj)))^(1/2))^2)*(T+Sij)/(T+Si);
% NANO PROPERTIES
mwn = 101.96     % 26.98(kg/kmol)
rhon = 3950 %2700;%(kg/m^3)% bulk
Cpn = 451 %   (10^3)*0.91;%(J/kgK)
volfM = pi()/6;% max volume fraction
Rnano =(1/2)*60e-9;% radius of avg nano particle in meters
df = 1.8;
Rg = (1/2)*250e-9;
ra =((Rg/Rnano)^df*Rnano^2)^(1/2);% from Keblinski Surface Area Equivilent radius
volf =@(T,Y) (rhon*(1/Y-1)/rhol(T)+1)^(-1);
inVolf = (Rg/Rnano)^(df-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4th page
MaxVolf = volfM*inVolf
%rhon = rhon*MaxVolf;% nano
rholMix =@(T,Yn) ((1-Yn)./rhol(T)+Yn./rhon).^-1;
%%
% ******************************************
% Simulation
% ******************************************
numdr = 500;
R = Rs(1);
dr = R./numdr;
r = 0:dr:Rs(1);
J = 0.0001*ones(loop,numdr);
RS = R/cos(pi/2-Theta);
h_max = -(sqrt(RS^2 - R^2)-RS);
h = RS*cos(asin(r(1:end-1)./RS)) - (RS-h_max);
Cap = (2*s(T0)/rhol(T0)/9.81)^0.5;
C = rhon*V0.*ones(1,size(r,2)-1);% Concentration from volume fraction
U = 0.*C;
% Volume representation
%V(1) = pi.*R^3.*Theta(1)./4;% Dunn et al. (Parabolic)
Beta = (1-cos(Theta))^2*(2+cos(Theta));
V(1) = (R/sin(Theta))^3*pi*(Beta)/3;% spherical from Erbil review 2012
Qtp = zeros(1,size(r,2)-1);
YpS = zeros(1,size(r,2)-1); % Mass Fraction of NANO
%Tls(1,:) = T0.*ones(1,size(r,2)-1);
count =0;
Yprev =YpS;
T = T0.*ones(1,size(r,2)-1);
dsdT =-(s(T0+.1) - s(T0))/(T0+.1 -T0);

if V0==0
    tf = V/7.37e-12;
else
    tf=V/5.4e-12;
end
M0=mean(C)*V;
for i=2:loop
    % SAVE VALUES
    if i == 2
        count=count+1;
        disp(['saved at ' num2str(i*dt) ' seconds'] );
        loc(count) = i;
        t(count) = 0;
        Tlsave(count,:) = T;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        Vsave(count) = V;
        V_dotsave(count) = 0;
        Thetasave(count) = Theta;
        Csave(count,:) = C;
        Usave(count,:) = U;
        
        VNsave(count) = sum(C.*2*pi.*r(1:end-1).*dr.*h);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5th page
        
        Cmodsave(count,:) = 1.*ones(1,size(r,2)-1);
    end
end
    
    % ****************************
    % Converge on a temperature
    % ****************************
    Tp = 0.99*T;
    Vp = V;
    Thetap = Theta;
    index_not_converged = ones(size(T,2),1);
    iloop = 0;
    converge = abs(T - Tp)./Tp;
    fprintf("\n")
    while converge(index_not_converged) > 0.0001 & iloop < 201
        fprintf("dummy ");
        V = Vp;
        Theta = Thetap;
        iloop = iloop+1;
        loc(count) = i;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 7th page
        
        t(count) = i*dt;
        Tlsave(count,:) = T;
        
        if iloop > 100
            loca = find(abs(T-Tp)./Tp > 0.01,1);
        end
        % Volume change
        %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        % VolFlux = 2*pi*r(1:end-1).*dr.*J./rhol(T);
        %V_dot = -2.*pi.*sum(J.*r(1:end-1).*dr./rhol(T));
        %%%%%%%%%%%%%%%%%%%%%%
        % Convert rate of change of volume to current volume
        % %%%%%%%%%%%%%%5555
%         V = V + V_dot.*dt;
%         if V < 0
%             disp(['Complete evap at ' num2str(i*dt) 'seconds of simulation time.'])
%             count=count+1;
%             disp(['saved at ' num2str(i*dt) ' seconds'] );
%             loc(count) = i;
%             %%%%%%%%%%%%%%%%%%%%%%55
%             t(count) = i*dt; 
%             Tlsave(count,:) = T;
%             Rsave(count,:) = R;
%             hsave(count,:) = h;
%             hh(count) = int_h;
%             Vsave(count) = V;
%             V_dotsave(count) = V_dot;
%             VolFlux_s(count,:) = VolFlux;
%             return
%         end
        % Evaluate Contact Angle, and height at new Volume
        % ******************************************
        int_theta = Theta;
        Theta = int_theta;
        int_h = -(sqrt(RS^2 - R^2)-RS);
        hprev = h;
        hh(count) = int_h;
        h = RS*cos(asin(r(1:end-1)./RS)) - (RS-int_h);
        delh = h - hprev;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        locmax = find(C > MaxVolf*rhon );
        if ~isempty(locmax)
            for invi =1:size(r,2)-1
                rloc = size(r,2)-invi;
                if C(rloc)>MaxVolf*rhon
                    if rloc-1 == 0
                        disp(['Max volume @ center @' num2str(MaxVolf),' Volume fraction, and ' num2str(i*dt) ' seconds' ])
                        V_dot = NaN;
                        V_return=V_dot;
                        return
                        
                    end
                end
            end
    end
end
figplot=1;
if figplot==1
    % Plotting
    % ******************************************
    loc = 1:count;
    figure(4)
    plot(t(loc),Thetasave(loc)*180/pi)
    hold on;
    plot(BF_t,BF_theta,'r');
    plot(NF_t,NF_theta,'g');
    grid on
    xlabel('Seconds');
    ylabel('Contact Angle (Degrees)')
    p=polyfit(t(loc),Thetasave(loc)*180/pi,1);
    plot(t(loc),polyval(p,t(loc)),'-.k','LineWidth',2);
    figure(9)
    plot(BF_t,BF_dtdt,'r');
    hold on
    plot(NF_t,NF_dtdt,'g');
    dtdt=gradient(180/pi*Thetasave(loc),t(loc));
    plot(t(loc),dtdt)
    p=polyfit(BF_t,BF_dtdt,1);
    plot(BF_t,polyval(p,BF_t),'-.k','LineWidth',2);
    xlabel('Seconds');
    ylabel('Contact Angle Rate (Degrees/sec)')
    h=figure(5)
    plot(t(loc),Vsave(loc),'LineWidth',2)
    hold on;
    plot(BF_t,BF_volume,'r')
    plot(NF_t,NF_volume,'g')
    pb=polyfit(BF_t,BF_volume,1);
    plot(BF_t,polyval(pb,BF_t),'r-.','LineWidth',2)
    pn=polyfit(NF_t,NF_volume,1);
    plot(NF_t,polyval(pn,NF_t),'g-.','LineWidth',2)
    grid on
    xlabel('Seconds');
    ylabel('Volume (m^3)')
    p=polyfit(t(loc),Vsave(loc),1);
    plot(t(loc),polyval(p,t(loc)),'-.k','LineWidth',2);
    legend(['Simulation, fit :' num2str(mean(V_dotsave(2:end)))],['Base Fluid, fit: ' num2str(pb(1)) ],['Nanofluid, fit: ' num2str(pn(1))]);
    print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_Evap_Rate.pdf'])
end