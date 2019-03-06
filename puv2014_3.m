clear;
close all
clc
%%
% This method is designed for shallow water area where the instruments
% could be exposed to the air; if not, the wave processing is the same, but
% alternal method of distinguishing different tides is needed.

 % hp height of the pressure sensor above the bottom, unit, m       
 % hv height of the velocity cell above the pressure sensor, unit, m 

time_interburst=300; % time interval between two bursts, unit: s
hp=0.646;
hv=-0.246;
dt=1/16;    %dt is the sample interval in s (typically 0.5 or 1 s
num=4096;numberneed=1*num;

numburst=numberneed/num;
% num: sample number of one burst
% numberneed: how much smaples used in compute one wave data, should be 2^N
% numberburst: how much bursts used in compute one wave data

FILE=load('CT_01_ADV5195.dat');          
burst=FILE(:,1);
p=FILE(:,15);  
vu=FILE(:,3);vv=FILE(:,4);vz=FILE(:,5);
SNR1=FILE(:,9);SNR2=FILE(:,10);SNR3=FILE(:,11);
LEN=length(burst); % number of total samples

cor1=FILE(:,12);cor2=FILE(:,13);cor3=FILE(:,14);

clear FILE
ini=burst(1);int=burst(num+1)-burst(1);ine=burst(LEN); 
index_burst=[ini:int:ine];
index_burst=index_burst(:);
LEN_burst=length(index_burst);
% ini: initial burst number
% int: interval number between two bursts
% ine: the last number of burst
% index_burst: the sequence of burst number
% LNE_burst: number of total bursts


%%
%--------------------------------------------------------------------------
% processing time
Time=load('CT_01_ADV5195.vhd');  
             
sec=Time(ini:int:ine,6);
min=Time(ini:int:ine,5);
hor=Time(ini:int:ine,4);
day=Time(ini:int:ine,2);
mon=Time(ini:int:ine,1);
yer=Time(ini:int:ine,3);
time=[yer mon day hor min sec];

time=datenum(time);                    
ttall=datestr(time,'mmm.dd,yyyy HH:MM:SS');
ttall=cellstr(ttall);
clear Time sec min hor day mon yer
%--------------------------------------------------------------------------
%%
% build structure
% in this code, index ii was used for burst number
clear DATA HEAD
for ii=1:int:LEN_burst  %
    index=(num*(ii-1)+1):num*ii;
    DATA(ii)=struct('pressure',p(index),'u',vu(index),'v',vv(index),...
        'cor1',cor1(index),'cor2',cor2(index),'cor3',cor3(index),...
        'SNR1',SNR1(index),'SNR2',SNR2(index),'SNR3',SNR3(index));
    % DATA was used to store the pressure, velocity, and correlation of
    %each burst
    
    % here to determine wether this burst data is good 检测数据质量，看相关性
    a=DATA(ii).cor1>=70&DATA(ii).cor1<101;
    b=DATA(ii).cor2>=70&DATA(ii).cor2<101;
    c=DATA(ii).cor3>=70&DATA(ii).cor3<101;
    d=DATA(ii).SNR1>30;
    e=DATA(ii).SNR2>30;
    f=DATA(ii).SNR3>30;
    g=a&b&c;gg=d&e&f;
    ggg=g&gg;
    LEN_gooddata=length(find(ggg));
    
    if LEN_gooddata/num>=0.985 
        % this is a strict condition that 99% samples in this burst should
        % be good, then it will be called good burst.如果一个burst里面有99%
        % 的数据是好的，认为这个是个好的burst。条件可以适当放宽到97%
        HEAD(ii)=struct('burstnumber',index_burst(ii),'time',time(ii),...
        'timestr',ttall(ii),'quality','good','tidenumber',nan,...
        'goodsamples',LEN_gooddata/num);
    else
        HEAD(ii)=struct('burstnumber',index_burst(ii),'time',time(ii),...
        'timestr',ttall(ii),'quality','bad','tidenumber',nan,...
        'goodsamples',LEN_gooddata/num);
    end
    %HEAD has five parameters: burstnumber, time, quality, tidenumber, and
    %    goodsamples
end
clear g gg ggg a b c d e f
clear cor1 cor2 cor3 vu vv vz SNR1 SNR2 SNR3
%%
for i=1:LEN_burst
    Gratio(i)=HEAD(i).goodsamples;
    Spressure(i)=mean(DATA(i).pressure);
end
figure;
subplot(2,1,1);
plot(1:LEN_burst,Gratio,'o');
ylabel('Gratio');
subplot(2,1,2);
plot(1:LEN_burst,Spressure,'ro');
ylabel('pressure recorded (dbar)');
%% --------------------------------------------------------------------------
% find the good burst number for next step
clear index
k=1;
for ii=1:LEN_burst
if strcmp(HEAD(ii).quality,'good')
index(k)=ii;        % the number of good burst will be stored in the index
k=k+1;
end
end
clear k

GOODburstnumber=length(index);              
%--------------------------------------------------------------------------
% divided bursts into different tides
% NOTE: this method need be changed if the instrument wasn't exposed to 
% the air between two tidal cyclings.
%%
k=1; % k is used to calculate the number of tidal cyclings.
for i=1:GOODburstnumber-1
     if (index(i+1)-index(i))>1&&(index(i+1)-index(i))<10
     disp('there is an abnormal discontinuity point during a tidal cycling')
     disp('at burst')
     disp(HEAD(index(i)).burstnumber)
     continue
     end
    
    HEAD(index(i)).tidenumber=k;
    tide(i)=k;
    % tide: tide number for varible 'index',also used for test
    
    if (index(i+1)-index(i))>40      % 滞空时间长达40个burst（约200分钟）表示已经经历了一个落干的时期，
       k=k+1;                        % 可以开始下一个潮周期的计算   
    end       
end
HEAD(index(GOODburstnumber)).tidenumber=k;
tide(GOODburstnumber)=k;

%% METHOD 3: 
for k=1:max(tide)
aa=find(tide==k);
ii=index(aa);
% calculate air pressure before the tide
    j=1;   
    pbefore(k)=mean(DATA(ii(1)-j).pressure);
% calculate air pressure after the tide
    j=1;
    pafter(k)=mean(DATA(ii(length(aa))+j).pressure);
% calculate p0
p0(k)=0.5*(pbefore(k) + pafter(k));
end
%%     
for i=1:GOODburstnumber
    
    xx=1:num;xx=xx(:);
    pp=polyfit(xx,DATA(index(i)).pressure,1);     
    trend=polyval(pp,xx);
    DATA(index(i)).wavepressure=DATA(index(i)).pressure-trend;%去除水位的趋势项，获取波浪的振幅
    k=HEAD(index(i)).tidenumber;      
    HEAD(index(i)).airpressure=p0(k); % acturally, this is not air pressue in this site
    HEAD(index(i)).waterdepth=polyval(pp,mean(xx))+hp...
        -HEAD(index(i)).airpressure;%取xx的平均值，这样得到的拟合值的误差最小
    % wave pressure data is stored in DATA variable
    % water depth is assigned to be the one at the MIDDLE of the
    % burst,取每个burst的中间水深作为表征水深
    % stored in HEAD variable
    % They only exist in "good" burst   
    
end
%%
% as structure data is not an array, may not be easily plotted or saved, so
% change it into array
for i=1:GOODburstnumber
Waterdepth(i)=HEAD(index(i)).waterdepth;% water depth
Time(i)=HEAD(index(i)).time; % time
end

figure
plot(Time,Waterdepth,'ro')
datetick('x','mmm.dd, HH:MM','keepticks') 
ylabel('red, pressure recorded (dbar）')

%% 
clear Time_wave Hs Waterdepth_wave peakF peakDir peakSpread tide_wave
% we have water depth data one per minute,now is five minute,one burst one
% depth
clear vu vv HEAD_wave DATA_wave index_wave tide_wave GOODburstnumber_wave
GOODburstnumber_wave=0;
 i=1;
 m=0;
for k=1:max(tide)
    aa=find(tide==k);
    ii=index(aa);
    for j=1:length(ii)
        time_seq(j)=HEAD(ii(j)).time;
    end
    tbegan=datestr(HEAD(ii(1)).time,21);
    
    for kk=1:length(aa)
        m=m+1;
        Waterdepthmid=HEAD(ii(kk)).waterdepth;        
        Time_wave=HEAD(ii(kk)).time;
    
        HEAD_wave(m)=struct('tide',k,'midwaterdepth',Waterdepthmid,...
       'midtime',Time_wave);
   
        vp(1+(m-1)*num:m*num)=DATA(ii(kk)).wavepressure;
        vu(1+(m-1)*num:m*num)=DATA(ii(kk)).u;
        vv(1+(m-1)*num:m*num)=DATA(ii(kk)).v;
        
        vp=vp(1:4096);vp=vp(:);
        vu=vu(1:4096);vu=vu(:);
        vv=vv(1:4096);vv=vv(:);
        
        DATA_wave(m)=struct('wavepressure',vp,'u',vu,'v',vv);
    end
    
end
GOODburstnumber_wave=m;

%% --------------------------------------------------------------------------
% as structure data is not an array, may not be easily plotted or saved, so
% change it into array


for i=1:GOODburstnumber_wave
Waterdepth_wave(i)=HEAD_wave(i).midwaterdepth;% water depth
Time_wave(i)=HEAD_wave(i).midtime; % time
tide_wave(i)=HEAD_wave(i).tide;
end

figure
plot(Time_wave,Waterdepth_wave,'ro')
title('figure of water depth with burst information')
ylabel('water depth, m')
xlabel('burst number')
datetick('x','mmm.dd, HH:MM')

%--------------------------------------------------------------------------
%% 编译到这里了
% calculate wave parameters

lf=.033;      %Hz - low frequency cutoff
hf=0.4; %Hz - high frequency cutoff
maxfac=200;   %   - maximum value of factor scaling pressure to waves
minspec=1e-4;  %m^2/Hz - minimum spectral level for computing
              %         direction and spreading
Ndir=0;     %deg - direction offset (includes compass error and 
              %      misalignment of cable probe relative to case
              % the offset for the Aquadopp Profiler is 0

parms=[lf hf maxfac minspec Ndir];
NF=128;

close all
clear Su Sp Dir Spread F dF
for i=1:GOODburstnumber_wave

[Su(:,i),Sp(:,i),Dir(:,i),Spread(:,i),F,dF,AA] = wds2014(DATA_wave(i).u,DATA_wave(i).v,DATA_wave(i).wavepressure,dt,NF,HEAD_wave(i).midwaterdepth,hp,hv,parms);
[Hs(i),peakF(i),peakDir(i),peakSpread(i)] = hs(Su(:,i),Sp(:,i),Dir(:,i),Spread(:,i),F,dF);% Copyright (C) 2001, Lee Gordon, NortekUSA LLC
[Hm(i),Ha(i),Ht(i)]=hm(DATA_wave(i).wavepressure);
       
end
figure
loglog(F,Sp)
  
figure
plot(Time_wave,peakF,'ro')

Time_wave=Time_wave(:);
Hs=Hs(:);
Waterdepth_wave=Waterdepth_wave(:);
peakF=peakF(:);
peakDir=peakDir(:);
peakSpread=peakSpread(:);
tide_wave=tide_wave(:);
save('RD2_data.mat','Time_wave','Hs','Waterdepth_wave','peakF','peakDir','peakSpread','tide_wave','Waterdepth_wave');
%     Hs           significant wave height
%     peakF        wave freauency (Hz) at the peak
%     peakDir      wave direction (deg) at the peak
%     peakSpread   wave spreading (deg) at the peak
xlswrite('RD2_0921-0925.xlsx',{'Time_wave','Hs','Waterdepth_wave','peakF','peakDir','peakSpread','tide_wave','Waterdepth_wave'},1,'a1');
xlswrite('RD2_0921-0925.xlsx',Time_wave,1,'a2');
xlswrite('RD2_0921-0925.xlsx',Hs,1,'b2');
xlswrite('RD2_0921-0925.xlsx',Waterdepth_wave,1,'c2');
xlswrite('RD2_0921-0925.xlsx',peakF,1,'d2');
xlswrite('RD2_0921-0925.xlsx',peakDir,1,'e2');
xlswrite('RD2_0921-0925.xlsx',peakSpread,1,'f2');
xlswrite('RD2_0921-0925.xlsx',tide_wave,1,'g2');
%%
figure
[AX,H1,H2]=plotyy(Time_wave,Hs,Time_wave,Waterdepth_wave,'plot');
set(get(AX(1),'Ylabel'),'String','Significant wave heights, m') 
set(AX(1),'XTickLabel',[])
set(AX(2),'XTickLabel',[])
set(get(AX(2),'Ylabel'),'String','Water depth, m') 
datetick('x','mmm.dd, HH:MM','keeplimits')
xlabel('Time') 
set(H1,'Marker','+','LineStyle', 'none')
set(H2,'Marker','o','LineStyle', 'none')

for k=1:max(tide_wave)
aa=find(tide_wave==k);
figure
[AX,H1,H2]=plotyy(Time_wave(aa),Hs(aa),Time_wave(aa),Waterdepth_wave(aa),'plot');
set(AX(1),'XTickLabel',[])
set(AX(2),'XTickLabel',[])
set(H1,'LineStyle', '-','Color','b')
set(H2,'LineStyle', '-','Color','k')
set(get(AX(1),'Ylabel'),'String','Significant wave heights, m') 
set(get(AX(2),'Ylabel'),'String','Water depth, m') 
datetick('x','mmm.dd, HH:MM','keeplimits')
xlabel('Time') 


end

%figure
%plotwds(Time_wave,Su,Sp,Dir,Spread,F,2);

figure
imagesc(Time_wave,log10(F),Dir,[-140 -80]);colorbar;
set(gca,'XTickLabel',[]);
title('Direction (deg N)');
ylabel('log(F) (Hz)');

figure
plot(Time_wave,peakDir)

figure
plot(Time_wave,peakF,'ro')

%%
% compare wave height with SBE
%Hs_SBE=WD(:,6);
%N_SBE=WD(:,11);

%Hs_SBE(N_SBE<150)=nan;

%figure
%plot(T_WD,Hs_SBE,Time_wave,Hs,'o')
%datetick('x','mmm.dd, HH:MM','keeplimits')



