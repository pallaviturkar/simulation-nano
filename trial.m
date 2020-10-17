BaseFluid = importdata('BF.xlsx', 1, 0);
BF_RADIUS = 0.0032/2;
BF_t = BaseFluid.data(4:end,1);% Seconds
loc = 50;
BF_t = BaseFluid.data(4:loc,1);% Seconds
BF_theta = BaseFluid.data(4:loc,2);%Degrees 
BF_dtdt = gradient(BF_theta,BF_t);
BF_Dradius = BaseFluid.data(4:loc,3)/1000;% mm to m
BF_height = BaseFluid.data(4:loc,4)/1000;% mm to m
BF_volume = BaseFluid.data(4:loc,5)/(1000^3);%   mm^3 to m^3 
BF_surfaceArea = BaseFluid.data(4:loc,6)/(1000^2);% mm^2 to m^2
BF_mass = BaseFluid.data(4:loc,7);% kg
BF_dmdt = gradient(BF_mass,BF_t);% kg/s
plot(BF_t,BF_volume);
xlabel('Time');
ylabel('Contact Angle')