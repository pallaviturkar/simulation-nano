T = 25+273;
Dg(T) = 2.38E-5;
Yinf =0;
Csat = 0.092;
R = 8.314;
J90 = Dg(T).*(1-Yinf).*Csat./R;% Units of Kg/m^2
lambda = 1/2 - Theta / pi;% found to be best match to FEM results
J0 = J90.*(0.27*Theta^2 + 1.30).*(0.6381 -0.2239*(Theta - pi/4)^2);
J = J0.*(1-r(1:end-1).^2./R.^2).^-lambda;
% loca = find(1./(1-r(1:end-1).^2./R.^2).^lambda > 7 ,1);
% J(loca:end) = J(loca).*exp(7 - 1./(1-r(loca:end-1).^2./R.^2).^lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iloop > 100
loca = find(abs(T-Tp)./Tp > 0.01,1);
J(loca:end) = J(loca).*exp(1./(1-r(loca).^2./R.^2).^lambda - 1./(1-r(loca:end-1).^2./R.^2).^lambda);
%disp('iloop greater than 100')
%sum(converge(index_not_converged))
end
%