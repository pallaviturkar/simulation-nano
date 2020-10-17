clear all
close all
dbstop if error
tic
df =[1.8 1.8];%[1.0,1.05,1.1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3];%[1.3 1.3 1.8 2 3];%
Rnano =(1/2)*60e-9;%(1/2).*[30:30:90].*1e-9;%
Rg =(1/2)*250e-9;%(1/2).*[70:50:300].*1e-9;
Sublayer= {'PTFE_Mod'};
skip=0;
for lk=1:length(Sublayer)
for i=1:size(df,2)
for j=1:size(Rnano,2)
for k=1:size(Rg,2)
if i==1
if skip ~= 1
disp('Running Pure fluid...')
[Vdot(i,j,k) Volf(i,j,k),bf,nf]=sessile_evap_v6(df(i),Rnano(j),Rg(k),0,Sublayer{lk});
Volf(i,j,k)=0;
else
[Vdot(i,j,k) Volf(i,j,k) bf]=[-7.3e-12 0 ];
end
else
[Vdot(i,j,k) Volf(i,j,k),bf,nf]=sessile_evap_v6(df(i),Rnano(j),Rg(k),0.02,Sublayer{lk});
end
end
end
end
end
h=figure(100*j)
plot(Volf(:,1,1),squeeze(-Vdot(:,1,1)),'+-');
hold on
if i>2
plot(0,-bf,'ko','LineWidth',2);
hold on;
end
if i==size(df,2)
plot([0 max(Volf)],[-nf -nf],'ro-','LineWidth',2);
end
plot(0:0.1:max(Volf),nf*ones(size(0:0.1:max(Volf)),1),'r+-')
axis([0 1 0 (10)*10^(-12)]);
xlabel('Maximum Volume Fraction')
ylabel('Evaporation rate')
legend('Simulation','Base-Fluid Sefiane and Bennacer @ 25C','Nanofluid Sefiane and Bennacer @ 25C')
grid on;
disp(['Nominal Evap of: ', num2str(Vdot(end,1,1))]);
disp('Percent Decrease of :');
PD=100*abs((Vdot(:,:,:)-Vdot(1,1,1))/Vdot(1,1,1))
print(h,'-dpdf',['absolute_study.pdf']);
h=figure(200)
plot(squeeze(Volf(:,:,:)),100-PD,'+-');
hold on;
plot(0:0.1:max(Volf),(nf/bf)*100*ones(size(0:0.1:max(Volf))),'r+-');
grid on;
xlabel('Max Volume Fraction');
ylabel('Percent of nominal evap')
legend('Simulation','Sefiane and Bennacer @ 25C')
print(h,'-dpdf',['percent_study.pdf']);
toc