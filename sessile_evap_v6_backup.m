function [V_return, MaxVolf, base_fluid,nano_fluid]=sessile_evap_v6(df, Rnano, Rg, V0, stringSublayer)
close all
clear all
%%
% ******************************************
% Load In Data from Sefiane and Bennacer
% @ 25 Celcius
% ******************************************
BaseFluid = importdata('Nanofluid.xlsx', 1, 0);
BF_RADIUS = 0.0012;
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
NanoFluid = importdata('NanoFluid.xlsx');
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
%%
% %%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% System
P =101325; % Pa
Yinf =0;%Gas Mass Fraction @ infinity
T0 = 273.15 + 25;%; 21.85 65
V0 = 0.02;
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
T0=298;  %%%%%My input
% Solid Surface parameters
% if strcmpi(stringSublayer,'Al')
% lams=237;
% elseif strcmpi(stringSublayer,'Steel')
% lams=50;
% elseif strcmpi(stringSublayer,'PTFE_Mod')
lams = 1;%.91;%0.25;%0.25;%237 ;% 237; % conductivity of
% aluminum from in W./mK
% elseif strcmpi(stringSublayer,'PTFE');
% lams =0.25;
% end
%0.91 matches well to base fluid results
hs = 0.0001;%0.0015;%1.5 mm thickness
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
mwn = 26.98    %101.96     % 26.98(kg/kmol)
rhon = 2700% 3950 %2700;%(kg/m^3)% bulk
Cpn =(10^3)*0.91    %451 %   ;%(J/kgK)
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
for i= 2:loop
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
        % ******************************************
        % Evaporation Mass Flux, J (Dunn et al 2007)
        % ******************************************
        Volf = C./rhon;% m^3 Al per m^3 solution
        % |kg nal | |m3 nAl | |m3 nal |
        % |-------|*|-------|=|-------|
        % |m3 tot | |kg nAl | |m3 tot |
        Y = Volf.*rhon./((1-Volf).*rhol(T) +Volf.*rhon);% kg al per kg solution
        Y = C./rholMix(T,Y);
        tY = C.*rhon./((1-C).*rhol(T) + Y.*rhon);
        Xlf = ((1-Y)./mwl)./(((1-Y)./mwl)+(Y./mwn));
        % Convert from mass fraction to mole
        %Xlf = ((1-YpS)./mwl)./(((1-YpS)./mwl)+(YpS./mwn));
        % Clausius
        Xgf =Xlf.*(101325/P).*exp((L(T,1)./(8.3144./mwl.*1000)).*(1./Tbf - 1./T));
        % Clausius Clayperion, gas phase mole frac of evap
        % Vapor Pressure 
        tXgf = Xlf.*(Pv(T)/P);
        Xgg = ((1-Yprev)./mwg)./((1-Yprev)./mwg+Yprev./mwl); % Gas phase mole fraction of nonevap
        if imag(Xgf) >0
            disp('Xgf is imag...');
            keyboard
        end
        Ygf = Xgf.*mwl./(Xgf.*mwl+Xgg.*mwg);
        % Mass fraction of fuel in vapor
        tYgf = Xgf.*mwl./(Xgf.*mwl+(1-Xgf).*mwg);
        % Mass fraction in vapor
        tCsat = (rhogf(T).*(1./Ygf-1)./transpose(rhog(T))+1).^(-1); % Concentration in vapor
        t2Csat =Ygf.*(1./rhogMix(T,Ygf));
        t3Csat = (1./(Xgf.*mwl + Xgg.*mwg)).*rhogMix(T,Ygf);
        Volfg = (rhogf(T).*(1./Ygf-1)./rhog(T)'+1).^(-1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%6th page
        
        
        %Csat =Volfg.*rhogf(T);
        % |m3 eth | |kg eth | |kg eth |
        % |-------|*|-------|=|-------|
        % |m3 tot | |m3 eth | |m3 tot |7
        Csat =Ygf.*(rhogMix(T,Ygf));% units of kg/m^3
        % |kg eth | |kg tot | |kg eth|
        % |-------|*|-------|=|------|
        % |kg tot | |m3 tot | |m3 tot|
        %Csat = 0.186 + 9.47*10^-3*(T-Tinf); % METHANOL VALUES
        %DgM =@(T) 1.5*10^(-5); % METHANOL VALUES
        %L =@(T,Y) 1155; % METHANOL VALUES
        %MASS FLUX FOR VARYING DEGREE CONTACT ANGLE
        J90 = Dg(T).*(1-Yinf).*Csat./R;% Units of Kg/m^2
        lambda = 1/2 - Theta / pi;% found to be best match to FEM results
        J0 = J90.*(0.27*Theta^2 + 1.30).*(0.6381 -0.2239*(Theta - pi/4)^2);
        J = J0.*(1-r(1:end-1).^2./R.^2).^-lambda;
        % loca = find(1./(1-r(1:end-1).^2./R.^2).^lambda > 7 ,1);
        % J(loca:end) = J(loca).*exp(7 - 1./(1-r(loca:end-1).^2./R.^2).^lambda);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iloop > 100
            loca = find(abs(T-Tp)./Tp > 0.01,1);
            J(loca:end) = J(loca).*exp(1./(1-r(loca).^2./R.^2).^lambda - 1./(1-r(loca:end-1).^2./R.^2).^lambda);
            %disp('iloop greater than 100')
            %sum(converge(index_not_converged))
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % Volume change
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VolFlux = 2*pi*r(1:end-1).*dr.*J./rhol(T);
        V_dot = -2.*pi.*sum(J.*r(1:end-1).*dr./rhol(T));
        %%%%%%%%%gf
        
        %%%%%%%%%%%%%
        % Convert rate of change of volume to current volume
        % %%%%%%%%%%%%%%5555
        V = V + V_dot.*dt;
        if V < 0
            disp(['Complete evap at ' num2str(i*dt) 'seconds of simulation time.'])
            count=count+1;
            disp(['saved at ' num2str(i*dt) ' seconds'] );
            loc(count) = i;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 7th page
            
            t(count) = i*dt;
            Tlsave(count,:) = T;
            Rsave(count,:) = R;
            hsave(count,:) = h;
            hh(count) = int_h;
            Vsave(count) = V;
            V_dotsave(count) = V_dot;
            VolFlux_s(count,:) = VolFlux;
            Jsave(count,:) = J;
            Csatsave(count,:) = Csat;
            Thetasave(count) = Theta;
            MaranSave(count) = Maran;
            Csave(count,:) =C;
            Usave(count,:) =U;
            return
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        int_theta = Theta;
        conver = .5;
        prev_theta = int_theta*.5;
        prev_vol = V*.5;
        while conver > 0.0001;
            RS = R/cos(pi/2-int_theta);
            int_h = -(sqrt(RS^2 - R^2)-RS);
            int_Vol = (pi*int_h/6)*(3*R^2+int_h^2);
            int_theta_n = int_theta - (int_Vol - V)*(int_theta-prev_theta)/(int_Vol-prev_vol);
            prev_vol = int_Vol;
            prev_theta = int_theta;
            int_theta = int_theta_n;
            conver = abs(int_Vol-V)/V;
        end
        %******************************************
        % Evaluate Contact Angle, and height at new Volume
        % ******************************************
        Theta = int_theta;
        hprev = h;
        h = RS*cos(asin(r(1:end-1)./RS)) - (RS-int_h);
        delh = h - hprev;
        %*******************************************
        % Esimate the local concentration as a result of Evap Flux J
        %*******************************************
        C = C.*hprev./h;
        %*******************************************
        % Estimate Velocity through change in volume
        %*******************************************
        %dVh = (delh./dt);%.*2*pi*r(1:end-1).*dr;
        %dVJ = J./rhol(T);%.*2*pi*r(1:end-1).*dr;
        %U = dVJ - dVh;
        %keyboard
        for rloc = 2:size(r,2)-1
            U(rloc) = (-1/rhol(T(rloc))./r(rloc)./h(rloc)).*sum(r(1:rloc).*dr.*(J(1:rloc) +rhol(T(1:rloc)).*delh(1:rloc)./dt));
        end
        tt=i*dt/tf; 
        rt=r(1:end-1)./R;
        u=(1/4).*(1./(1-tt)).*1./rt.*((1-rt.^2).^-lambda - (1-rt.^2));
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
        Cprev=C;
        for rloc = 2:size(r,2)-1
            C(rloc) = (C(rloc-1)*U(rloc-1)*dt*2*pi*r(rloc-1)*h(rloc-1)- C(rloc) *U(rloc) *dt*2*pi*r(rloc) *h(rloc)+ C(rloc) *dr *2*pi*r(rloc) *h(rloc))/(2*pi*r(rloc)*h(rloc)*dr);
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
        totError =M0 - sum(C.*pi.*r(1:end-1).*2*dr.*h)...
            /(sum(2*pi.*r(1:end-1).*h.*dr))*V;
        C(end) =C(end) + totError/(2*pi*r(end-1)*dr*h(end-1));
        % ******************************************
        % Marangoni Calc
        % ******************************************
        Maran = mean(dsdT.*L(T,1).*Dg(T).*(1-Yinf).*Csat./Theta./s(T)./lamff(T));
        % ******************************************
        % Gravity Calc
        % ******************************************
        Grav = (9.81*(rhon-rhol(T))/rhon)./(18/4*visfl(T).*(delh./dt)./ra^2);
        % ******************************************
        % Spatially varying Temp (Dunn et al. 2007)
        % ******************************************
        % Added Heat conduction/diffusive transport in Qair
        % Bm = (Ygf-Yinf)./(1-Ygf);
        Aff =A(T,visf(T),visf(T),mwl,mwl,Sf,Sf,Sff);
        Agg =A(T,visg(T),visg(T),mwg,mwg,Sg,Sg,Sgg);
        Agf =A(T,visg(T),visf(T),mwg,mwl,Sg,Sf,Sgf);
        Afg =A(T,visf(T),visg(T),mwl,mwg,Sf,Sg,Sfg);
        % lamMix = Xgf.*lamf(T)./(Xgf.*Aff+Xgg.*Afg)+Xgg.*lamg(T)./(Xgf.*Agf+Xgg.*Agg);
        % rhoMix = ((Ygf)./rhogf(T)+(1-Ygf)./(rhog(T)')).^1;
        % Le =((rhoMix.*CpMix(T,Ygf).*Dg(T)./lamMix)).^-1;
        % CpgBAR = (Cpg(T)+Cpg(Tinf))./2;% air avg. specific heat
        % CpfBAR =(2.*Cpf(T)+Cpf(Tinf))./3;%avg. gas phase fuel specific heat
        % phi = (CpfBAR./CpgBAR).*(1./Le);
        % Bt = ((1+Bm).^phi)-1;
        % Qair = rhol(T).*V_dot.*(CpfBAR.*(Tinf-T)./Bt);
        %
        % Tn = Tinf - ((L(T,1).*J + Qair).*(h./lamff(T) + hs./lams));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Tgrad = T0.*(1-r(1:end-1) /R*.01);
        %Tinf = Tgrad;
        Tn = Tinf - ((L(T,1).*J).*( h./lamff(T) + hs./lams));
        loca = find(r(:) > R/2 ,1);
        e=0:(h(loca))/100:h(loca);
        Tfuns =@(z) Tinf - L(T(loca),1)*J(loca)/lams*(z+hs);
        Tfun =@(z) Tinf -L(T(loca),1)*J(loca)*(z./lamff(T(loca)) + hs./lams);
        BL(1,:) = e;
        BL(2,:) = Tfun(e);
        BL(3,:) = -hs:hs/100:0;
        BL(4,:) = Tfuns(BL(3,:));
        Tinf_forced = ((L(T,1).*J).*( h./lamff(T) + hs./lams)) + 298;
        loc=find(isnan(Tn));
        loc1=find(Tn<0);
        if ~isempty(loc) | ~isempty(loc1) | ~isempty(find(Tn>Tinf))
            disp('Temperature acting all funky...')
            beep
            keyboard;
            loc=find(Tn>Tinf | isnan(Tn) | Tn<0)
            Tn(loc)=Tn(loc-1);
        end
        Ts = Tinf - L(T,1).*J/lams*(hs);
        Qsurface = lams*(Ts-Tinf)/hs;
        %**********************************************************
        % Neglect unsteady temperature via hu et al 2005 Appendix A
        % *********************************************************
        % look at Stanton Number/
        Tp = T;
        T = Tn;
        clear index_not_converged
        index_not_converged = find(abs(T - Tp)./Tp > 0.0001);
        converge = abs(T - Tp)./Tp;
    end
    % ******************************************
    % PREVIOUS VALUES
    % ******************************************
    Yprev = Ygf;
    % ******************************************
    % SAVE VALUES
    % ******************************************
    %if mod((dt*i),interval) ==0 | i==10
        count=count+1;
        disp(['saved at ' num2str(i*dt) ' seconds'] );
        disp([' mean (Csat) ' num2str(mean(Csat))]);
        disp([' mean (tCsat) ' num2str(mean(tCsat))]);
        disp([' mean (t2Csat) ' num2str(mean(t2Csat))]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([' mean (t3Csat) ' num2str(mean(t3Csat))]);
        disp([' mean (Ygf) ' num2str(mean(Ygf))]);
        disp([' mean (tYgf) ' num2str(mean(tYgf))]);
        disp([' mean (Xgf) ' num2str(mean(Xgf))]);
        disp([' mean (tXgf) ' num2str(mean(tXgf))]);
        disp([' mean (Xlf) ' num2str(mean(Xlf))]);
        disp([' Partial pressure = Xgf (using Pv)' num2str(mean(Pv(T)./P))]);
        loc(count) = i;
        t(count) = i*dt;
        Tlsave(count,:) = T;
        Tinf_forcesave(count,:) =Tinf_forced;
        Rsave(count,:) = R;
        hsave(count,:) = h;
        hh(count) = int_h;
        Vsave(count) = V;
        V_dotsave(count) = V_dot;
        VolFlux_s(count,:) = VolFlux;
        Jsave(count,:) = J;
        Csatsave(count,:) = Csat;
        Ygfsave(count,:) = Ygf;
        Thetasave(count) = Theta;
        MaranSave(count) = Maran;
        GravSave(count,:) = Grav;
        Csave(count,:) = C;
        Xlfsave(count,:) = Xlf;
        Ysave(count,:) = Y;
        Usave(count,:) = U;
        U1save(count,:) = U1;
        VNsave(count) = sum(C.*2*pi.*r(1:end-1).*dr.*h);
        Cmodsave(count,:) = Cmod;
        BL_s(count,:,:) = BL;
   % end
end
%disp(['Vdot of ' num2str(V_dot) ' @ ' num2str(numdr) ' Radial Steps, and dt of ' num2str(dt)', 'Seconds'])
figplot=1;
if figplot==1
    % ******************************************
    % Plotting
    % ******************************************
%     loc = 1:count;
%     figure(4)
%     plot(t(loc),Thetasave(loc)*180/pi)
%     hold on;
%     plot(BF_t,BF_theta,'r');
%     plot(NF_t,NF_theta,'g');
%     grid on
%     xlabel('Seconds');
%     ylabel('Contact Angle (Degrees)')
%     p=polyfit(t(loc),Thetasave(loc)*180/pi,1);
%     plot(t(loc),polyval(p,t(loc)),'-.k','LineWidth',2);
%     figure(9)
%     plot(BF_t,BF_dtdt,'r');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hold on
%     plot(NF_t,NF_dtdt,'g');
%     dtdt=gradient(180/pi*Thetasave(loc),t(loc));
%     plot(t(loc),dtdt)
%     p=polyfit(BF_t,BF_dtdt,1);
%     plot(BF_t,polyval(p,BF_t),'-.k','LineWidth',2);
%     xlabel('Seconds');
%     ylabel('Contact Angle Rate (Degrees/sec)')
%     h=figure(5)
%     plot(t(loc),Vsave(loc),'LineWidth',2)
%     hold on;
%     plot(BF_t,BF_volume,'r')
%     plot(NF_t,NF_volume,'g')
%     pb=polyfit(BF_t,BF_volume,1);
%     plot(BF_t,polyval(pb,BF_t),'r-.','LineWidth',2)
%     pn=polyfit(NF_t,NF_volume,1);
%     plot(NF_t,polyval(pn,NF_t),'g-.','LineWidth',2)
%     grid on
%     xlabel('Seconds');
%     ylabel('Volume (m^3)')
%     p=polyfit(t(loc),Vsave(loc),1);
%     plot(t(loc),polyval(p,t(loc)),'-.k','LineWidth',2);
%     legend(['Simulation, fit :' num2str(mean(V_dotsave(2:end)))],['Base Fluid, fit: ' num2str(pb(1)) ],['Nanofluid, fit: ' num2str(pn(1))]);
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_Evap_Rate.pdf']);
%     figure(6)
%     plot(t(loc),abs(V_dotsave(loc))/(1e-12),'o')
%     hold on;
%     p=polyfit(t(loc),abs(V_dotsave(loc))/(1e-12),1);
%     %plot(BF_t,polyval(p,BF_t),'-.','LineWidth',2);
%     plot(BF_t,-BF_dmdt/790/(1e-12),'r');
%     plot(NF_t,-NF_dmdt/790/(1e-12),'g');
%     p=polyfit(BF_t,-BF_dmdt/790/(1e-12),1);
%     plot(polyval(p,BF_t),'-.r','LineWidth',2);
%     p=polyfit(NF_t,-NF_dmdt/790/(1e-12),1);
%     plot(polyval(p,NF_t),'-.g','LineWidth',2);
%     %plot(-gradient(polyval(pv,BF_t))/(1e-12),'g-.','LineWidth',2)
%     grid on
%     grid
%     xlabel('Seconds');
%     ylabel('- Volume rate of change (nl/s)')
%     grid
%     disp(['Vdot of ' num2str(V_dot) ' @ ' num2str(numdr) ' Radial Steps, and dt of ' num2str(dt) ' Seconds']);
%     disp([' Vdot mean of ' num2str(mean(V_dotsave(2:end)))]);
%     figure(7)
%     plot(R,abs(V_dotsave(loc))/(1e-12),'o')
%     hold on;
%     plot(R,-BF_dmdt/790/(1e-12),'r+')
%     plot(R,-NF_dmdt/790/(1e-12),'g+')
%     grid on
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     grid
%     xlabel('Radius (m)');
%     ylabel('- Volume rate of change (nl/s)')
%     axis([0 2e-3 0 30 ])
%     grid
%     figure(10)
%     plot(R,(abs(dtdt(2))./(abs(V_dotsave(loc(2)))/(1e-12))),'o');
%     grid on
%     hold on
%     xlabel('Initial radius')
%     ylabel('Rate of change of theta over rate of change of volume')
%     figure(11)
%     plot(R,(abs(dtdt(2))),'o');
%     grid on
%     hold on
%     plot(R,BF_dtdt,'r+')
%     xlabel('Initial radius')
%     ylabel('Rate of change of theta')
%     figure(8)
%     plot(t(loc),MaranSave(loc))
%     hold on
%     plot([t(loc(1)) t(loc(end))],[0.05 0.05],'r');
%     axis([0 i*dt 0 1])
%     grid
%     xlabel('Seconds');
%     ylabel('Instantanious Marangoni number')
%     
%     tic
%     figure(20)
%     plot(t(loc),VNsave(loc))
%     hold on;
%     plot(0,M0,'ro');
%     grid on
    for j =1:count
        c='';
        % if j==1
        % c='r';
        % elseif j==size(loc,2)
        % c='k';
        % end
        if V0==0
            c='';
        else
            c='r';
        end
        figure(1)
        plot(r(1:end-1),hsave(j,:),c)
        hold on
        figure(2)
        plot(r(1:end-1),Tlsave(j,:)-273.15,c)
        hold on
        %%%%%%%%%%%%%%%%%%%%%%
        figure(3)
        plot(r(1:end-1)./R,Jsave(j,:)*1000,c)
        hold on;
        figure(12)
        plot(r(1:end-1)./R,VolFlux_s(j,:),c);
        hold on;
        figure(13)
        plot(r(1:end-1)./R,Csave(j,:)/rhon,c);
        hold on;
        figure(17)
        plot(r(1:end-1)./R,Ysave(j,:),c);
        hold on;
        % figure(18)
        % plot(r(1:end-1)./R,GravSave(j,:),c);
        % hold on;
        % figure(14)
        % plot(r(1:end-1)./R,Usave(j,:),c);
        % hold on;
        %
        figure(15)
        plot(r(1:end-1)./R,Csatsave(j,:),c);
        hold on;
        figure(19)
        plot(r(1:end-1)./R,Csave(j,:),c);
        hold on;
        plot(r(1:end-1)./R,Ygfsave(j,:),c);
        hold on;
        figure(21)
        plot(r(1:end-1)./R,Tinf_forcesave(j,:),c);
        hold on;
        grid on
        figure(16)
        plot(r(1:end-1)./R,Usave(j,:),c)
        hold on;
        plot(r(1:end-1)./R,U1save(j,:),[''])
%         if j>1
%             figure(22)
%             plot(squeeze(BL_s(loc(j),2,:))- 273.15,squeeze(BL_s(loc(j),1,:)),'LineWidth',2);
%             hold on
%             plot(squeeze(BL_s(loc(j),4,:))- 273.15,squeeze(BL_s(loc(j),3,:)),'LineWidth',2);
%             grid on
%             hold on
%             xlabel('Temperature');
%             ylabel(['Height @ Half Radius and time' num2str(t(loc(j)))]);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         end
    end
%     h=figure(1)
%     xlabel('Radius (m)');
%     ylabel('Height (m)');
%     title(['Evaporation every ' num2str(interval) ' seconds'])
%     str(1)={['Fractal Dimension: ' num2str(df)]};
%     str(2)={['Nanoparticle Radius: ' num2str(Rnano) ' m']};
%     str(3)={['Agglomerate Radius: ' num2str(Rg) ' m']};
%     str(4)={['Initial volume fraction: ' num2str(V0)]};
%     str(5)={['Maximum Volume fraction: ' num2str(MaxVolf)]};
%     text(0.1e-3,0.1*hsave(1,1),str,'BackgroundColor','w');
%     grid on
%     axis equal
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_height.pdf']);
%     h=figure(2)
%     xlabel('Radius (m)');
%     ylabel('Temperature (Celcius)');
%     title(['Evaporation every ' num2str(interval) ' seconds'])
%     text(.1e-3,Tlsave(end,1)-273.15+2,str,'BackgroundColor','w');
%     grid on
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_Temp.pdf']);
%     h=figure(3)
%     xlabel('Normalized Radius (n.d.)');
%     ylabel('Evaporative Flux J (g/m^2s)');
%     title(['Evaporation every ' num2str(interval) ' seconds'])
%     text(0.1,3,str,'BackgroundColor','w');
%     grid on
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_Flux.pdf']);
%     figure(12)
%     xlabel('Normalized Radius (n.d.)');
%     ylabel('m^3/s Evaporation Volume per time at each radial station')
%     grid on
%     h= figure(13)
%     hold on
%     plot(r(1:end-1)/R,ones(size(r(1:end-1)/R))*MaxVolf,'r-.');
%     xlabel('Normalized Radius (n.d.)');
%     ylabel('Local Volume Fraction');
%     text(0.1,.16,str,'BackgroundColor','w');
%     grid on
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_VolFrac.pdf']);
%     h=figure(17)
%     xlabel('Normalized Radius (n.d.)');
%     ylabel('Local Mass Fraction')
%     text(0.1,.16,str,'BackgroundColor','w');
%     grid on
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_MassFrac.pdf']);
%     h=figure(16)
%     title(['Evaporation every ' num2str(interval) ' seconds'])
%     xlabel('Normalized Radius (n.d.)');
%     ylabel('Local Radial Velocity')
%     text(0.1,0.8e-5,str,'BackgroundColor','w');
%     grid on
%     print(h,'-dpdf',['fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_Vel.pdf']);
%     figure(15)
%     xlabel('Normalized Radius (n.d.)');
%     ylabel('Saturation Concentration')
%     grid on
%     figure(19)
%     xlabel('Normalized Radius (n.d.)');
%     ylabel(' Gas Phase Fuel mass fraction at the surface (Ygf)');
%     grid on
end
% V_dotsave
% base_fluid= pb(1);
% nano_fluid = pn(1);
% V_return = mean(V_dotsave(2:end));
% %V_return = V_dot;
% p_test=polyfit(t,Vsave,1);
% disp([' Slope fit of : ' num2str(p_test(1))]);
% saveas(gcf, 'figure 1.png') %'figure 2.png','figure 3.png','figure 4.png','figure 5.png','figure 6.png','figure 7.png','figure 8.png','figure 9.png','figure 17.png','figure 22.png');
% %save(['sublayer_' stringSublayer '_fractal_dimension' num2str(df) '_VolF' num2str(V0) num2str(Rg) '_MassFrac.mat']);
return
end
