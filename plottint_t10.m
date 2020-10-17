% close all
% clear all
% %%%%%%%%%%%%%
% NanoFluid = importdata('NanoFluid.xlsx', 1, 0);
%Theta
clear all
close all
clc
tic

%% Directory
dir_main = 'C:\Users\Pallavi\Desktop\Nano_modelling';
cd(dir_main)
dir_data=fullfile(dir_main,'Data');

%%

dir_path=fullfile(dir_data,'comparison_NF.xlsx');
Data = xlsread(dir_path, 'Sheet1', 'A2:AJ16500'); %%%% reads the specified worksheet.
% figure(1)
%  NF_t =Data(:,3);
% NF_theta= Data(:,4);
% NF_t1 = Data(:,11);
% NF_theta_model =Data(:,12);
%  plot(NF_t,NF_theta,'r','LineWidth',2);
%  hold on
% plot(NF_t1,NF_theta_model,'k','LineWidth',2);
% xlabel('time (sec)');
% ylabel('Contact angle (degree)');
% legend('Paper data', 'modeling data');
% %Volume
% figure(2)
% NF_t2 = Data(:,5);
% NF_vp= Data(:,6);
% NF_t1 = Data(:,11);
% NF_v_model = Data(:,13);
% plot(NF_t,NF_vp,'r','LineWidth',2);
% hold on
% plot(NF_t1,NF_v_model,'k','LineWidth',2);
% xlabel('time (sec)');
% ylabel('Volume');
% legend('Paper data', 'modeling data');

% % Evaporation flux, J
% figure(3)    
 NF_r1 = Data(:,1);
 NF_Jp= Data(:,2);
 NF_r = Data(:,3);
NF_J_code =Data(:,4);
 plot(NF_r1,NF_Jp,'r','LineWidth',2);
 hold on
 plot(NF_r,NF_J_code,'k','LineWidth',2);
xlabel('Normalized Radius');
ylabel('Evaporative Flux J (g/m^2s)');
legend('Paper data-t10', 'modeling data-t10');

%  Drop height
%  figure(4)    
%  NF_r2 = Data(:,20);
%  NF_hp= Data(:,21);
%  NF_h_code = Data(:,23);
%  plot(NF_r2,NF_hp,'r','LineWidth',2);
%  hold on
%  plot(NF_r,NF_h_code,'k','LineWidth',2);
%  xlabel('Normalized Radius');
%  ylabel('Drop height (mm)');
%  legend('Paper data', 'modeling data');

%  %  Drop Temperature
%  figure(5)    
%  NF_r3 = Data(:,6);
%  NF_Tp= Data(:,7);
% % NF_rcode = NanoFluid.data(:,8);
%  NF_T_code =Data(:,9);
%  plot(NF_r3,NF_Tp,'r','LineWidth',2);
%  hold on
%  plot(NF_r,NF_T_code,'k','LineWidth',2);
%  xlabel('Normalized Radius');
%  ylabel('Drop temperature (K)');
%  legend('Paper data', 'modeling data');
 
%   %  Solid Fraction
  figure(6)    
 NF_r4 = Data(:,11);
 NF_Volf_P= Data(:,12);
 %NF_r1code = Data(:,35);
 NF_Volf_code = Data(:,14);
 plot(NF_r4,NF_Volf_P,'r','LineWidth',2);
 hold on
 plot(NF_r,NF_Volf_code,'k','LineWidth',2);
 xlabel('Normalized Radius');
 ylabel('Solid fraction');
 legend('Paper data', 'modeling data');
    
% %  RAdial velocity
%  figure(7)    
%  NF_r5 = Data(:,29);
%  NF_U_P= Data(:,30);
%  %NF_r1code = NanoFluid.data(:,13);
%  NF_U_code = Data(:,31);
%  plot(NF_r5,NF_U_P,'r','LineWidth',2);
%  hold on
%  plot(NF_r,NF_U_code,'k','LineWidth',2);
%  xlabel('Normalized Radius');
%  ylabel('Solid fraction');
%  legend('Paper data', 'modeling data');