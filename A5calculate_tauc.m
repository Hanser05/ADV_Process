% Calculate shear stress
% 用分离波浪之后的turbulence
clear all;clc;close all
data = load('adv_turbulence.mat');
burst = data.newdata(:,1); time = data.newdata(:,2);
d = 0.2+0.246; % distance from pressure sensor to seabed
Depth = data.newdata(:,3) + d; 
u_turbulence = data.newdata(:,4); 
v_turbulence = data.newdata(:,5); 
w_turbulence = data.newdata(:,6);
density = 1018; % sea water density
xmax= datenum([2017 12 7 12 0 0]);  % end time
xmin= datenum([2017 12 1 0 0 0]);  % start time
xstep = datenum([0 0 1 0 0 0]);  
%
data2 = load('advClean_despiking_majorminor.mat');
u_major = data2.new_data1(:,4); 
v_minor = data2.new_data1(:,5);
w = data2.new_data1(:,6);
U = sqrt(u_major.^2+v_minor.^2);
%
gap = find(diff(burst)~=0);
gap = [0;gap;length(burst)];
for i = 1:length(gap)-1
    index=gap(i)+1:gap(i+1);    
    DATA(i)=struct('burst',burst(index),'time',time(index),'depth',Depth(index),'U',U(index),...
        'u_major',u_major(index),'v_minor',v_minor(index),'w',w(index),...
        'u_turbulence',u_turbulence(index),'v_turbulence',v_turbulence(index),...
        'w_turbulence',w_turbulence(index),'uw_turbulence',u_turbulence(index).*w_turbulence(index),...
        'u_tst',u_major(index)-mean(u_major(index)),'v_tst',v_minor(index)-mean(v_minor(index)),...
        'w_tst',w(index)-mean(w(index)));
end
clear i
%% shear stress
for i = 1:length(gap)-1
   tauc_rs(i,1) = -density*mean(DATA(i).uw_turbulence);  % RS
   tauc_tke(i,1) = (mean(DATA(i).u_turbulence.^2)+mean(DATA(i).v_turbulence.^2)+...
       mean(DATA(i).w_turbulence.^2))/2*0.19*density;  %TKE
   
   tauc_tke_tst(i,1) = (mean(DATA(i).u_tst.^2)+mean(DATA(i).v_tst.^2)+...
       mean(DATA(i).w_tst.^2))/2*0.19*density;  %TKE
   
   tauc_mtke(i,1) = density*0.9*mean(DATA(i).w_turbulence.^2); %m-TKE
   ave_time(i,1) = mean(DATA(i).time);
   ave_depth(i,1) = mean(DATA(i).depth);
   ave_u_major(i,1) = mean(DATA(i).u_major);
   ave_v_minor(i,1) = mean(DATA(i).v_minor);
   % ID
%    RawW = DATA(i).w_turbulence;
   RawW = DATA(i).w;
   N = length(RawW);
   Fs = 8;
  [Spectral,F]=pwelch(RawW-mean(RawW),[],[],N,Fs,'oneside');%Welch
   u_star(i,1) = (2*pi*0.4*0.2/mean(DATA(i).U))^(1/3)*(max(F.^(5/3).*Spectral)/0.68)^0.5; 
%    u_star1(i,1) = (2*pi*0.4*0.25/mean(s_tur))^(1/3) * (trapz(f,p)/FS/0.69*2)^0.5;% 16/2=8, 

   tauc_id(i,1) = density*u_star(i,1)^2;
end
ave_U = sqrt(ave_u_major.^2+ave_v_minor.^2);
tauc_id(tauc_id(:,1)>5,:)= nan;
%%
interval = datenum([0 0 0 0 30 0]);
blank= find(diff(ave_time)>interval);
blank = [0;blank;length(ave_time)];
figure;
gca1 = subplot('position',[0.1 0.1 0.8 0.35]);%% 
tauc_rs= abs(tauc_rs);
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    % DATA = struct('t',time_tauc(index),'id',tauc_id(index));
    plot(ave_time(index),tauc_rs(index),'k.-','linewidth',2,'markersize',10); hold on;
    plot(ave_time(index),tauc_tke(index),'r.-','linewidth',2,'markersize',10);hold on;
plot(ave_time(index),tauc_mtke(index),'g.-','linewidth',2);hold on;
plot(ave_time(index),tauc_id(index),'b.-','linewidth',2);hold on;
end
xmax= datenum([2017 5 7 12 0 0]);
xmin= datenum([2017 5 10 0 0 0]);
xstep = datenum([0 0 1 0 0 0]);  
%

fz = 15;
set(gca,'fontsize',fz,'fontname','times new roman','linewidth',1.5);
xlabel('Time (day/hour)');
ylabel('\tauc (N/m^2)');
legend('RS','TKE','TKEw','ID');legend('boxoff');
box off;
set(gca1,'fontsize',fz,'fontname','times new roman','linewidth',1.5);
%set(gca1,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]); 
datetick('keepticks');
%set(gca1,'ylim',[0 2],'xtick',[0:1:2]); 
%set(gca,'PlotBoxAspectRatio',[6 4 1]); 
 parameters = [ave_time tauc_rs tauc_tke tauc_mtke tauc_id];
xlswrite('τc_CT01.xlsx',{'ave_time' 'tauc_rs' 'tauc_tke' 'tauc_mtke' 'tauc_id'},'sheet1','a1');
xlswrite('τc_CT01.xlsx',parameters,'sheet1','a2');




