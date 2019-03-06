%% Step6: calculate wave parameters
tic; clear all; clc;
data = load('advClean_despiking.mat');
load('advClean_wavepressure.mat');
burst = data.new_data2(:,1);
time = data.new_data2(:,2); % time
vp = data.new_data2(:,3); % raw pressure
vu = data.new_data2(:,4); % eastern velocity
vv = data.new_data2(:,5); % northern velocity
depth = data.new_data2(:,11);
%
xmax= datenum([2017 5 10 12 0 0]);  % start time
xmin= datenum([2017 5 17 0 0 0]);  % end time
xstep = datenum([0 0 1 0 0 0]);  
gap = find(diff(burst)~=0); gap = [0;gap;length(burst)];
for i = 1:length(gap)-1
    index = gap(i)+1:gap(i+1);
    DATA(i) = struct('time',time(index),'burst',burst(index),'pressure',vp(index),'u',vu(index),'v',vv(index),...
        'depth',depth(index));
end
%% PUV and sta to cal wave parameters
%
lf = .025;      %Hz - low frequency cutoff
hf = 1; %Hz - high frequency cutoff
hf2 = 2;  % Hz - resolution frequece
maxfac = 200;   %   - maximum value of factor scaling pressure to waves
factor = 0.4;  %   
minspec = 1e-4;  %m^2/Hz - minimum spectral level for computing
              %         direction and spreading
Ndir=0;       %deg - direction offset (includes compass error and 
              %      misalignment of cable probe relative to case
              % the offset for the Aquadopp Profiler is 0
NF=128;
parms=[lf hf maxfac minspec Ndir];
hp = 0.246+0.2; % hp height of the pressure sensor above the bottom, unit, m  
hv = -0.396; % hv height of the velocity cell above the pressure sensor, unit, m
dt = 1/16;  
lambda = 10000;
%
j = 1;
for i = 1:length(gap)-1
    w0 = mean(DATA(i).depth);
    Time(j,1) = mean(DATA(i).time);
    Vu(j,1) = mean(DATA(i).u);
    Vv(j,1) = mean(DATA(i).v);
     w_depth = w0; 
    p_temp = DATA(i).pressure; 
    p = kick_trend(p_temp,lambda);

    [Hm,Ha,H3,wT,HT,mT,H10,L]=sta_wave(p,dt,w_depth,hp,hf2,factor);
    % PUV
    rd_u = DATA(i).u; rd_v = DATA(i).v;
    [Su,Sp,Dir,Spread,F,dF,DOF,AA] = wds2014(rd_u,rd_v,p,dt,NF,w_depth,hp,hv,parms);
    [Hs,peakF,peakDir,peakSpread] = hs(Su,Sp,Dir,Spread,F,dF);
    waterdepth(j,1) = w_depth;
    
    [Hm_a(i),Ha_a(i),Ht_a(i)]=hm(wavepressure(i).wavepressure);
    %
    Hm_sta(j,1) = Hm; 
    Tm_sta(j,1) = wT;
    Ha_sta(j,1) = Ha; 
    Ta_sta(j,1) = HT; 
    Hs_sta(j,1) = H3; 
    Ts_sta(j,1) = mT;
    L_sta(j,1) = L;   
    H10_sta(j,1) = H10;
    peakDir_puv(j,1) = peakDir; 
    peakF_puv(j,1) = peakF;
    Hs_puv(j,1) = Hs;
     
    
    j = j+1;
end
%
figure
plot(Time,Hs_puv,'r');
hold on
plot(Time,Hs_sta,'g');
% plot(Time,Hs_spe,'b');
legend('PUV','Statistic');
% cal wave number
K_sta = wavek(1./Tm_sta,waterdepth); %refer to Whitehouse 2000, T=Tp, Tp is the peak period
                                            % w = 2*pi*f = 2*pi/T = 2*pi/Tp
                                            % so f is the peak frequency
K_puv = wavek(peakF_puv,waterdepth);
%% Bottom velocity from simple linear wave theroy
Uw_puv = pi*Hs_puv/sqrt(2).*peakF_puv./sinh(waterdepth.*K_puv);
Uw_sta = pi*Hs_sta/sqrt(2)./Tm_sta./sinh(waterdepth.*K_sta);
%% cal Uw(Bottom velocity from near-bottom velocity measurments)
%(Wiberg and Sherwood 2008)
% spectra analyze
Fs = 16;
for i = 1:length(gap)-1
    RawU = DATA(i).u;
    N_u = length(RawU);
    [S_RawU,f]=pwelch(RawU-mean(RawU),[],[],N_u,Fs,'oneside');
    turU = PhaseDec(RawU,Fs);
    [S_turU,f] = pwelch(turU-mean(turU),[],[],N_u,Fs,'oneside');
    S_waveU = S_RawU-S_turU;
    %
    RawV = DATA(i).v;
    N_v = length(RawV);
    [S_RawV,f]=pwelch(RawV-mean(RawV),[],[],N_v,Fs,'oneside');
    turV = PhaseDec(RawV,Fs);
    [S_turV,f] = pwelch(turV-mean(turV),[],[],N_u,Fs,'oneside');
    S_waveV = S_RawV-S_turV;
    %
    iw = find(f>=0.1 & f<=0.25); %wave frequency
    Spectral = S_waveU + S_waveV;
    Uw_spe(i,1) = sqrt(abs(2*sum(Spectral(iw)*0.002))); 
    Hs_spe(i,1) = 4*sqrt(sum(Spectral*0.002/(2*pi*cosh(0.3*K_puv(i,1))/sinh(K_puv(i,1)*waterdepth(i,1)))^2));
end

%%

interval = datenum([0 0 0 0 30 0]);
blank= find(diff(Time)>interval);
blank = [0;blank;length(Time)];
 figure;
 subplot('position',[0.1 0.6 0.8 0.35]);
 plot(data.new_data2(:,2),depth+0.4,'k.-','linewidth',1,'markersize',10);
 datetick; 
 xlabel('time(day)');
 ylabel('Depth');
 title('Depth');
 box off;

gca1 = subplot('position',[0.1 0.1 0.8 0.35]);  %

for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    % DATA = struct('t',time_tauc(index),'id',tauc_id(index));
    plot(Time(index),Uw_puv(index),'k.-','linewidth',2,'markersize',10); hold on;
    plot(Time(index),Uw_sta(index),'r.-','linewidth',2,'markersize',10);hold on;
plot(Time(index),Uw_spe(index),'g.-','linewidth',2);
end

datetick; 
%set(gca1,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]); 
legend('Uw-PUV','Uw-sta','Uw-spectra');legend('boxoff');
set(gca,'linewi',1.5,'fontsize',14,'fontname','times new roman');
%set(gca,'PlotBoxAspectRatio',[6 4 1]); 
xlabel('time(day)');
ylabel('Wave orbital velocity (m/s)');
title('GT01-Uw');
box off;
%%
Time_str = cellstr(datestr(Time));
wavedir = peakDir_puv;
for i = 1:length(wavedir)
    if wavedir(i) < 0
        wavedir(i) = wavedir(i)+360;
    end
end
Hm_a= Hm_a';
Ha_a= Ha_a';
Ht_a= Ht_a';
parameters = [Time waterdepth Hs_sta Hs_puv Tm_sta 1./peakF_puv wavedir K_sta K_puv Uw_sta Uw_puv Uw_spe];
parameters2 =[Hm_a Ha_a Ht_a ];
parameters3 =[Time Hm_sta Hs_sta Ha_sta H10_sta];
xlswrite('adv_wave_parameters.xlsx',...
    {'time','depth','Hs_Sta','Hs_puv','Tp_sta','Tp_puv','wavedir','k_sta','k_puv','Uw_sta','Uw_puv','Uw_spe'},'sheet1','a1');
xlswrite('adv_wave_parameters.xlsx',...
    parameters,'sheet1','a2');
xlswrite('adv_wave_parameters.xlsx',...
    {'Hm_a','Ha_a','Ht_a'},'sheet1','ag1');
xlswrite('adv_wave_parameters.xlsx',...
    parameters2,'sheet1','ag2');
xlswrite('adv_wave_parameters.xlsx',...
    {'time','Hm_sta', 'Hs_sta' ,'Ha_sta' ,'H10_sta'},'sheet2','a1');
xlswrite('adv_wave_parameters.xlsx',...
    parameters3,'sheet2','a2');
toc;

















