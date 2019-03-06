% Step7: cal tauw
clear all;
clc;
d50 = 12.074*10^-6;  
%Hs=data(:,8);Tz=data(:,7);height=data(:,2);
rd_wave = xlsread('d2_CT01_wave_parameters.xlsx');
Time= rd_wave(:,1);Hs=rd_wave(:,4);Tz=rd_wave(:,6);height=rd_wave(:,2);
density = 1027;  
viscosity = 1.353830*10^-6;  
%for i=1:length(data)
for i=1:length(rd_wave)
    [touWave,Uwave]=touWaveProcedure(Hs(i),Tz(i),height(i),d50,density,viscosity);
    tou(i,1)=touWave;
    u(i,1)=Uwave;
end
tauC= xlsread('¦Óc_adv.xlsx');
tauc_rs= tauC(:,2);
tauc_tke= tauC(:,3);
tauc_mtke= tauC(:,4);
tauc_id= tauC(:,5);
data = load('advClean_despiking.mat');
velocity_dir= data.aveUVW(:,4);
%% 
tau = [tou tauc_rs tauc_tke tauc_mtke tauc_id velocity_dir];
xlswrite('adv_wave_parameters.xlsx',{'touWave','tauc_rs','tauc_tke','tauc_mtke','tauc_id','velocity_dir'},'sheet1','N1');
xlswrite('adv_wave_parameters.xlsx',tau,'sheet1','N2');