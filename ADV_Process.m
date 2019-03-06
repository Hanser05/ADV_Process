clear;clc;
load('adv_enu.mat');
ADVdata=adv_data;clear adv_data;
%{
load('adv4841_transformed.dat');
ADVdata(:,3:5)=adv4841_transformed; clear adv4841_transformed;
%}
load('N3181212.vhd');
ADV5049_vhd=N3181212;


%%
Year=ADV5049_vhd(:,3);
Mon=ADV5049_vhd(:,1);
Day=ADV5049_vhd(:,2);
Hour=ADV5049_vhd(:,4);
Min=ADV5049_vhd(:,5);
Sec=ADV5049_vhd(:,6);
Time_vhd=datenum([Year Mon Day Hour Min Sec]);
Time_vhd=[Time_vhd,ADV5049_vhd(:,7)];
%%
N=8192; %numbers
fs=16; %frequency /Hz
hv=1.146;%distance from presssure probe to seabed
dis=hv-0.246-0.15;%distance from velocity probe to seabed
gap=find(diff(ADVdata(:,1))>0.5);
gap=[0;gap;size(ADVdata,1)];
%process of PSD
%cycle
%ADV5049_Tur=zeros(size(ADVdata,1),5)*nan;
Tc_Rey=zeros(length(gap)-1,1)*nan;%Reynolds shear stress
Tc_TKE=zeros(length(gap)-1,1)*nan;
Tc_TKEw=zeros(length(gap)-1,1)*nan;
Ustar_ID=zeros(length(gap)-1,1)*nan;%Reynolds shear stress
Uc_dir=zeros(length(gap)-1,1)*nan;%direction of main current for each burst
TKEc=zeros(length(gap)-1,1)*nan;%TKE of turbulence
TKEcw=zeros(length(gap)-1,1)*nan;%TKE of fluctuation items
Uw=zeros(length(gap)-1,1)*nan;%wave orbital velocity
%Urms=zeros(length(gap)-1,1)*nan;%RMS of Turbulence
Tc_With_Uw=zeros(length(gap)-1,1)*nan;
Umean=zeros(length(gap)-1,1)*nan;%velocity at east
Vmean=zeros(length(gap)-1,1)*nan;%wave at north
Wmean=zeros(length(gap)-1,1)*nan;%wave at verticle
Dep=zeros(length(gap)-1,1)*nan;%water depth
Burst=zeros(length(gap)-1,1)*nan;%Bursts
Time=zeros(length(gap)-1,1)*nan;%Time
tol=10^-6;%values for iter
for i=1:length(gap)-1;
    TempData=ADVdata(gap(i)+1:gap(i+1),:);    
    SNR1=TempData(:,9);
    SNR2=TempData(:,10);
    SNR3=TempData(:,11);
    Cor1=TempData(:,12);
    Cor2=TempData(:,13);
    Cor3=TempData(:,14);
    a = SNR1>=20;
    b = SNR2>=20;
    c = SNR3>=20;
    d = Cor1>=80;
    e = Cor2>=80;
    f = Cor3>=80;
    g=a&b&c&d&e&f;
    LEN_gooddata=length(find(g));
    if LEN_gooddata/N>0.8;

        tempU=TempData(:,3);%velocity at east
        tempV=TempData(:,4);%velocity at north
        tempW=TempData(:,5);%velocity at verticle
        %despiking  u                   
                    U0 = tempU;
                    [tempU,spike] = advPrepare(tempU);    
                    while abs(U0-tempU) >= tol;
                        U0 = tempU;
                        tempU = advPrepare(tempU);
                    end                   
                    %despiking  v    
                    V0 = tempV;
                   [tempV,spike] = advPrepare(tempV);                
                   while abs(V0-tempV) >= tol;
                       V0 = tempV;
                       tempV = advPrepare(tempV);
                   end
                   %despiking  w1
                   W0 = tempW;
                   [tempW,spike] = advPrepare(tempW);
                   while abs(W0-tempW) >= tol;
                       W0 = tempW;
                       tempW = advPrepare(tempW);
                   end

        mU=mean(tempU);%mean velocity at east
        mV=mean(tempV);%mean velocity at north
        mW=mean(tempW);%mean velocity at verticle
        uvw=sqrt(mU^2+mV^2+mW^2);
        dir=atan2(mU,mV);%angle of main stream direction
        UL=tempU*sin(dir)+tempV*cos(dir);%velocity at main stream direction
        VL=-tempU*cos(dir)+tempV*sin(dir);%velocity at secondary stream direction
        turU=tempU-mU;%turbulence at east
        turV=tempV-mV;%turbulence at north
        turW=tempW-mW;%turbulence at verticle
        turUL=UL-mean(UL);%turbulence at horizontal main stream 
        turVL=VL-mean(VL);
        turU_dwave=PhaseDec(turU,fs);
        turV_dwave=PhaseDec(turV,fs);
        turUL_dwave=PhaseDec(turUL,fs);%delete wave orbital velocity
        turVL_dwave=PhaseDec(turVL,fs);%delete wave orbital velocity
        turW_dwave=PhaseDec(turW,fs);%
        %UW=turU_dwave.*turW;% u'w' to delete wave orbital
        %VW=turV_dwave.*turW;% v'w' to delete wave orbital
        
        %ADV5049_Tur(gap(i)+1:gap(i+1),:)=[TempData(:,1),TempData(:,2),U0,V0,W0];
%Reynolds
        Tc_Rey(i)=-1023*mean(turUL_dwave.*turW_dwave);
        Uc_dir(i)=dir;%main stream direction
%---TKE shear stress------
        Tc_TKE(i) = (mean(turUL_dwave.^2)+mean(turVL_dwave.^2)+mean(turW_dwave.^2))/2*0.19*1023; %TKE   
%---TKEw shear stress------
        Tc_TKEw(i) = 0.9*1023*mean((turW_dwave).^2); %TKEw 

%ID method
%       erate(i)=(1/a*Sww_f_mean2).^(3/2)*2*pi/abs(Ubar(i));
        Fl = uvw/dis; % lower limit 
        Fu = fs/2;  % upper limit    
        
        % turb1d
        [f1,s1]=turb1d_ef(turW1_dwave,uvw1,fs,4,N/4,2);
        F1=f1(f1<Fu&f1>Fl);S1=s1(f1<Fu&f1>Fl);
        tmp1=F1.^(5/3).*S1;
        Ustar_IDw1(j,i)=(2*pi*0.4*dis1/uvw1)^(1/3)*(mean(tmp1)/0.69)^0.5;
        Touc_ID_smooth(j,i)=density*Ustar_IDw1(j,i).^2;
        epsin_smooth(i)=(nanmean(tmp1)/0.69)^1.5*(2*pi/uvw);
        
        % pwelch
        [s2,f2]=pwelch(turW1_dwave,[],[],N,fs,'oneside');
        F2=f2(f2<Fu&f2>Fl);S2=s2(f2<Fu&f2>Fl);
        tmp2=F2.^(5/3).*S2;
        Ustar_IDw2(j,i)=(2*pi*0.4*dis1/uvw1)^(1/3)*(mean(tmp2)/0.69)^0.5;
        Touc_ID_pw(j,i)=density*Ustar_IDw2(j,i).^2;
        epsin_pw(i)=(nanmean(tmp2)/0.69)^1.5*(2*pi/uvw);
        %[s2,f2]=pwelch(turW2_dwave,[],[],N,fs,'oneside');
        %[s1,f1]=pwelch(turW,[],[],N,16,'oneside');
        %F=f1(f1<=3&f1>=1.5);
        %S=s1(f1<=3&f1>=1.5);%Lissa J. MacVean, frequency:1.5-3Hz
        %tmp = F.^(5/3).*S;   
        %Tc_ID(i) = (2*pi*0.4*dis/uvw)^(1/3)*(mean(tmp)/0.678)^0.5;
%Tc do not delete waves        
        uw=mean(turU.*turW);vw=mean(turV.*turW);
        Tc_With_Uw(i)=1023*sqrt(uw^2+vw^2);
%TKE
        TKEc(i)=mean(turU_dwave.^2)+mean(turV_dwave.^2)+mean(turW.^2);
        TKEcw(i)=mean(turU.^2)+mean(turV.^2)+mean(turW.^2);
        

        Uw(i,:)=ADV_Uw(tempU,tempV,1/15,1,fs);
        %Urms(i)=sqrt(mean(turU.^2+turV.^2));
        Dep(i)=mean(TempData(:,15))+hv;
        Umean(i)=mU;
        Vmean(i)=mV;
        Wmean(i)=mW;
        %Burst(i)=mean(TempData(:,1));
        %Time(i)=Time_vhd(Time_vhd(:,2)==Burst(i),1);   
    end
    Burst(i)=mean(TempData(:,1));
    Time(i)=Time_vhd(Time_vhd(:,2)==Burst(i),1); 
    id1=[];id2=[];
end

ADV_Analyze=struct('Tc_Rey',Tc_Rey,...
    'Uc_dir',Uc_dir,...
    'Tc_TKE',Tc_TKE,...
    'Tc_TKEw',Tc_TKEw,...
    'Ustar_IDw1',Ustar_IDw1,...
    'Ustar_IDw2',Ustar_ID,...
    'Touc_ID_pw',Touc_ID_pw,...
    'Touc_ID_smooth',Touc_ID_smooth,...
    'epsin_pw',epsin_pw,...
    'epsin_smooth',epsin_smooth,...
    'Tc_With_Uw',Tc_With_Uw,...
    'TKEc',TKEc,...
    'TKEcw',TKEcw,...
    'Uw',Uw,...
    'Umean',Umean,...
    'Vmean',Vmean,...
    'Wmean',Wmean,...
    'Dep',Dep,...
    'Time',Time,...
    'Epsin',epsin,...
    'Burst',Burst);
save ('ADV_Analyze_v2.mat','ADV_Analyze');                                                                                                                                                                                                                                                                          

figure;
plot(ADV_Analyze.Time,ADV_Analyze.Tc_Rey,'r');
hold on;
plot(ADV_Analyze.Time,ADV_Analyze.Tc_TKE,'b');
plot(ADV_Analyze.Time,ADV_Analyze.Tc_TKEw,'g');
plot(ADV_Analyze.Time,ADV_Analyze.Ustar_ID.^2*1023,'k');
datetick;
legend('Reynolds','TKE','TKEw','ID');
legend('Reynolds','TKE','TKEw','ID');