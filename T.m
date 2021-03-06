function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
%clear all
% ******************************************
% Load In experimental data
% ******************************************
NanoFluid = importdata('NanoFluid.xlsx', 1, 0);
NF_RADIUS = 0.0024/2;
loc = 55;
NF_t = NanoFluid.data(1:loc,1);                  % Seconds
NF_theta = NanoFluid.data(1:loc,2);
%dtdt = gradient(BF_theta,t);
% height = BF.data(1:loc,4)/1000;            % mm to m
% volume = BF.data(1:loc,5)/(1000^3);        % mm^3 to m^3
% surfaceArea = BF.data(4:loc,6)/(1000^2);   % mm^2 to m^2

%%%%%%%%%%%%%%%%%
% System
P =101325;                              % Pa
Yinf =0;                                % Gas Mass Fraction @ infinity
T0 = 273.15 + 25;%; 21.85 65
Ru =8.314;
tf=40;
Tinf = 273.15 + 25;
V0= 0.02;
lams =1;
tss = 0.0015;
Rs(1) = NF_RADIUS;


%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.001;
seconds = floor(NF_t(end));
loop = seconds./dt;
interval = 10; %

%%%%%%%%%%%%
% Liquid Parameters
Dg = 1.24E-05;
rhol = 789;
L = 9.2040 * 100000;   %821.92*1000;
lamff = 0.171;
mwl = 46;              % Molecular weight of the liquid
Tb = 78.8+273;         % Boiling point of the fuel;
rhogf =0.1453;
%%%%%%%%%%%%%%%%
%Gas properties
mwg =28.97;                   % kg/kmol
rhog = 1.18;

% NANO PROPERTIES
mwn = 26.58             %101(kg/kmol)
rhon = 2700             %3950;%(kg/m^3)% bulk
Cpn = (10^3)*0.91       %451 %   (10^3)*0.91;%(J/kgK)
volfM = pi()/6;         % max volume fraction
Rnano =(1/2)*60e-9;     % radius of avg nano particle in meters
df = 1.7;
Rg = (1/2)*250e-9;
% ra =((Rg/Rnano)^df*Rnano^2)^(1/2);% from Keblinski Surface Area Equivilent radius
inVolf = (Rg/Rnano)^(df-3);
MaxVolf = volfM*inVolf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4th page

%rhon = rhon*MaxVolf;       % nano
%rholMix =@(T,Yn) ((1-Yn)./rhol(T)+Yn./rhon).^-1;

% Simulation
numdr = 500;
R = Rs(1);
dr = R./numdr;
temp = 0:dr:Rs(1);
r = temp(1:end-1);
count =0;

% Ts = T0.*ones(1,size(r,2)-1);

% <START>======Initial values==========
j = 1;   % j is time parameter
Theta(1) = NF_theta(1)*pi/180; % radians, NF_theta is coming from excel sheet
Beta = (1-cos(Theta)).^2.*(2+cos(Theta));
% J in Kg/m^2/s
%J(j, :) = 0.0001*ones(j,size(r,2));
% Volume 
V(j) =  (R/sin(Theta)).^3*pi*(Beta)/3;
% Concentration in Kg/m^3
C(j, :) = rhon*V0.*ones(1,size(r,2));
% Volume fraction
volf(j, :) = C(j, :)./rhon;   %m^3 Al per m^3 solution
% Liquid velocity in meter
U(j, :) = 0.*C;
% Radius of spherical cap
RS(j) = R/cos(pi/2-Theta);
% Droplet height in meter
h_max = -(sqrt(RS.^2 - R.^2)-RS);
h(j, :) = RS*cos(asin(r(1:end)./RS)) - (RS-h_max);
% Drop surface temperature in Kelvin
Tls(j, :) = T0.*ones(1,size(r,2));
% initial Mass fraction of nano-particle
YpS = zeros(1,size(r,2));
Yprev =YpS;
% saturated vapour concentration
Csat(j, size(r)) = 0;
M0=mean(C)*V;
% <END>======Initial values==========


for j =  2:loop   % time loop
    
    % Calculation of saturated vapour concentration
    Y = volf(j-1,:).*rhon./((1-volf(j-1, :)).*rhol + volf(j-1, :).*rhon);
    Xgg = ((1-Yprev)./mwg)./((1-Yprev)./mwg+Yprev./mwl);    % Gas phase mole fraction of nonevap
    Xlf = ((1-Y)./mwl)./(((1-Y)./mwl)+(Y./mwn));
    Xgf =Xlf.*(101325/P).*exp((L./(8.3144./mwl.*1000)).*(1./Tb - 1./Tls(j-1, :)));
    Ygf = Xgf.*mwl./(Xgf.*mwl+Xgg.*mwg);
    rhogMix =Ygf.*rhogf +(1-Ygf).*rhog';                    % linear mixture rule
    Csat(j, :) = Ygf.*(rhogMix);                % units of kg/m^3
    
    
    % Evaporation flux
    lambda(j) = 1/2 - Theta(j-1)/pi;
    J0 = (Dg.*( 1-Yinf).*Csat(j, :)./R).*(0.27*Theta(j-1).^2 + 1.30).*( 0.6381 -0.2239.*(Theta(j-1) - pi/4).^2);
    J(j, :) = J0.*(1-r(1:end).^2./R.^2).^-lambda(j);
    
    % volume flux Vf
    VolFlux(j, :) = 2*pi*r(1:end).*dr.*J(j,:)./rhol;
    V_dot(j, :) = -2.*pi.*sum(J(j, :).*r(1:end).*dr./rhol);
    V(j) =  V(j-1) + V_dot(j, :).*dt;
    
    if V < 0  % Break the time loop once we get negative volume and go for ploting
        t(count) = j*dt;
        break
    end
    
    % Calculation of theta
    c1 = (3*V(j))/(pi*R^3);
    syms x
    eqn = (((1-cos(x))^2)*(2+cos(x)))/(sin(x))^3 == c1;
    S = solve(eqn, x, 'Real', true);
    Theta(j) = double(S);
    
    % Calculation of h at new volume and theta
    RS(j) = R/cos((pi/2)-Theta(j));
    disp("RS " +RS(j));
    int_h = -(sqrt(RS(j)^2 - R^2)-RS(j));
    % hprev = h(j-1,:);
    h(j, :) = RS(j)*cos(asin(r(1:end)./RS(j))) - (RS(j)-int_h);
    delh = h(j, :) - h(j-1, :);
    
    % Temperature
    Tls(j, :) = Tinf - ((L.*J(j, :)).*( h(j, :)./lamff + tss./lams));
    Ts  = Tinf- ((L.*J).*( (h+tss)./lams));   % Ignore
    
    % Esimate the local concentration as a result of Evap Flux J
        %*******************************************
        C(j,:) = C(j-1,:).*h(j-1,:)./h(j,:);
        %*******************************************
        % Estimate Velocity through change in volume
        %*******************************************
        %dVh = (delh./dt);%.*2*pi*r(1:end-1).*dr;
        %dVJ = J./rhol(T);%.*2*pi*r(1:end-1).*dr;
        %U = dVJ - dVh;
        %keyboard
        for rloc = 2:size(r,2)
            U(j, rloc) = (-1/rhol)./r(rloc)./h(j,rloc).*sum(r(1:rloc).*dr.*(J(j, 1:rloc) +rhol.*delh(1:rloc)./dt));
        end
        tt=j*dt/tf; 
        rt=r(1:end-1)./R;
        u=(1/4).*(1./(1-tt)).*1./rt.*((1-rt.^2).^-lambda(j) - (1-rt.^2));
        U1=u.*R./tf;
        U=U1;
        U(1) = 0;
        loc=find(isnan(U));
        if ~isempty(loc)
            keyboard
            U(loc)=0;
        end
        % ******************************************
        % Modify Concentration due to convection
        % ******************************************
       for rloc = 2:size(r,2)
        C(j,rloc) = (dt*2*pi*(C(j-1,rloc-1)*U(j, rloc-1)*r(rloc-1)*h(j-1,rloc-1)-...
           C(j-1,rloc) * U(j, rloc) *r(rloc) *h(j-1,rloc))+...
            dr *2*pi*C(j-1,rloc)*r(rloc) *h(j-1,rloc))/(2*pi*r(rloc)*h(j,rloc)*dr);
        if rloc == 2
            C(j, 1) = C(j, 2);  
        end
       end
        Cmod = C./Cprev;
        if ~isempty(find(C<0))
            keyboard
            
            
        end
        % *******************************************************
        % Find if concentration is over max, then place inward
        %********************************************************
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
                    %disp('Max Volf...')
                    Excess = 2*pi*r(rloc)*dr*h(rloc)*(C(rloc)-MaxVolf*rhon);
                    C(rloc) = MaxVolf*rhon;
                    C(rloc-1) = C(rloc-1) + Excess/(2*pi*r(rloc-1)*dr*h(rloc-1));
                    totError=M0 - sum(C.*pi.*r(1:end-1).*2*dr.*h).../(sum(2*pi.*r(1:end-1).*h.*dr))*V;
                        %%%%%%%%%%%%%%%%%%%%%%%%
                    C(rloc-1) =C(rloc-1) + totError/(2*pi*r(rloc-1)*dr*h(rloc-1));
                    if C(rloc-1)<0
                        disp('C less than 0...')
                        beep
                        keyboard
                    end
                end
            end
        end
        totError =M0 - sum(C.*pi.*r(1:end).*2*dr.*h)...
            /(sum(2*pi.*r(1:end).*h.*dr))*V;
       C(j,end) =C(j,end) + totError(j)/(2*pi*r(end-1)*dr*h(j,end-1));
    %  Volume fraction
    volf(j, :) = C(j,:)./rhon;
    Yprev = Ygf;
    %%%%%%%%%%%%%%%%%%%%%%%%
    %if mod((dt*j),interval) ==0 | i==10
    count=count+1;
    disp(['saved at ' num2str(j*dt) ' seconds'] );
    disp([' mean (Csat) ' num2str(mean(Csat))]);
    loc(count) = j;
    % Time data for Volume and theta plot
    t(j) = j*dt;
    
end



% Plotting
loc = 1:j;

count = j;
for j =1:count
    c='';
    figure(3)
    plot(r(1:end-1)./R,h(j,:)*1000,c)
    xlabel('Normalized Radius');
    ylabel('Height (m)');
    hold on
    figure(4)
    plot(r(1:end-1)./R,J(j,:)*1000,c)
    xlabel('Normalized Radius');
    ylabel('Evaporative Flux J (g/m^2s)');
    hold on
    figure(5)
    plot(r(1:end-1)./R,u(j,:),c);
    xlabel('Normalized Radius');
    ylabel(' Radial velocity');
    hold on
    figure(6)
    plot(r(1:end-1)./R,Tls(j,:),c);
    hold on
    figure(7)
    plot(r(1:end-1)./R,Ts(j,:),c)
    hold on
    figure(8)
    plot(r(1:end-1)./R,C(j,:),c)
    hold on
    figure(9)
    plot(r(1:end-1)./R,volf(j,:),c)
    hold on
end
