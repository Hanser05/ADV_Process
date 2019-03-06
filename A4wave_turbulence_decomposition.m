%% Step4: wave-turbulence decomposition using PhaseDec function
% ESA(Energy spectrum analysis)
clear all; clc; close all;
data = load('advClean_despiking.mat');
burst = data.new_data1(:,1); time = data.new_data1(:,2);
d = 0.2+0.246; % distance from pressure sensor to seabed
Depth = data.new_data1(:,3) + d; 
u_major = data.new_data1(:,4); %uMajor
v_minor = data.new_data1(:,5); %vMinor
w = data.new_data1(:,6);
%% struct to save data
gap = find(diff(burst)~=0);
gap = [0;gap;length(u_major)];
for ii = 1:length(gap)-1
    index=gap(ii)+1:gap(ii+1);    % uu is u', uuww is u'w'
    DATA(ii)=struct('time',time(index),'burst',burst(index),'u_major',u_major(index),'v_minor',v_minor(index),...
        'w',w(index),'uu',u_major(index)-mean(u_major(index)),'depth',Depth(index));
end
%% plot decomposition 
figure
RawU1 = DATA(1).u_major;
N = length(RawU1);
Fs = 16;
[Spectral1,F1]=pwelch(RawU1,[],[],N,Fs,'oneside');%Welch 
loglog(F1,Spectral1,'b','LineWidth',2.5); hold on
%
% RawU2 = PhaseDec(RawU1,Fs);
% [Spectral2,F2]=pwelch(RawU2,[],[],N,Fs,'oneside');%Welch 
% loglog(F2,Spectral2,'r','LineWidth',2)
x = 10^(-0.5):10^(-1):7;
loglog(x,x.^(-5/3)*10^(-3.5),'k','LineWidth',2.5);
% legend('wave+turbulence spectrum','turbulence spectrum (ESA)','f^{-5/3}');
legend('boxoff')
fz = 20;
set(gca,'fontsize',fz,'fontname','times new roman','linewidth',1.5);
set(gca,'PlotBoxAspectRatio',[6 4 1]); %设置x轴和y轴的比为6:4
xlabel('Frequency (Hz)');
ylabel('S_u (m^2/s)');
%% plot decomposition
figure
RawU1 = DATA(50).u_major;
N = length(RawU1);
Fs = 16;
[Spectral1,F1]=pwelch(RawU1-mean(RawU1),[],[],N,Fs,'oneside');%Welch 
loglog(F1,Spectral1,'k','LineWidth',2.5); hold on
%
RawU2 = PhaseDec(RawU1,Fs);
[Spectral2,F2]=pwelch(RawU2-mean(RawU2),[],[],N,Fs,'oneside');%Welch 
loglog(F2,Spectral2,'r','LineWidth',2.5)
x = 10^(-0.5):10^(-1):7;
loglog(x,x.^(-5/3)*10^(-3.5),'k','LineWidth',2.5);
legend('wave+turbulence spectrum','turbulence spectrum (ESA)','f^{-5/3}');
legend('boxoff')
fz = 20;
set(gca,'fontsize',fz,'fontname','times new roman','linewidth',1.5);
set(gca,'PlotBoxAspectRatio',[6 4 1]); 
xlabel('Frequency (Hz)');
ylabel('S_u (m^2/s)');
%% wave-turbulence decomposition
for i = 1:length(gap)-1
    index1 = gap(i)+1:gap(i+1); 
    u_turbulence(index1) = PhaseDec(DATA(i).u_major-mean(DATA(i).u_major),Fs);
    v_turbulence(index1) = PhaseDec(DATA(i).v_minor-mean(DATA(i).v_minor),Fs);
    w_turbulence(index1) = PhaseDec(DATA(i).w-mean(DATA(i).w),Fs);
end 
%%  save
newdata = data.new_data1;
newdata(:,4) = u_turbulence;
newdata(:,5) = v_turbulence;
newdata(:,6) = w_turbulence;
save('adv_turbulence.mat','newdata');


