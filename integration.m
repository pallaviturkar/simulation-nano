clear all
close all
clc
tic
%% Directory
dir_main = 'C:\Users\Pallavi\Desktop\Nano_modelling';
cd(dir_main)

dir_data = dir_main; 
syms x y
% fun = x;
% 
% syms Js J(y) R
% int (int ((1-x/R)^2)^-(0.5 - y/pi),'x', 0, 3)
fnum = sprintf('BF.xlsx');
dir_path = fullfile(dir_data, fnum);
Data = xlsread(dir_path);
r =  Data(:,3);
Theta= Data(:,2);

% r=[3.131760079, 3.328019722 3.308354039, 3.209628432,3.170694255];
% Theta=[39.22479824, 35.81679399, 34.42544055, 34.90831623, 35.35739593];


D =2.38*1e-5; C = 1.92; Y = 0; R =3.221;
f = @(x, y)((D.*C.*(1-Y)).*(0.27.*y.^2 + 1.3).*(0.6381 - 0.2239.*(y - pi./4).^2)./R).*(1- (x./R).^2).^(-(1./2-y./pi));
p = []
J = []
for i = 1:length(r)
    J(i) =round( integral2(f, 0,r(i), 0, Theta(i)), 4);
    fprintf("J=%0.4f\n", J(i));   
    p(i)= r(i)/R;
end
 plot (p, -J)
 xlabel('\it r/R','fontsize',20);
ylabel('\it Evaporative flux','fontsize',20);



 