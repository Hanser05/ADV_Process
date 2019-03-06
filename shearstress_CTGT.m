clear all;
clc;
data1 = xlsread('d2_CT01_wave_parameters.xlsx');   
time1 = data1(:,1);            
tauW1= data1(:,14);
tauC1= data1(:,15);  %rs
tauCW_soulsby1= data1(:,22);
tauCW_grant1= data1(:,23);
xmin = datenum([2017 5 10 0 0 0]);  %-------------------------
xmax = datenum([2017 5 17 0 0 0]);  %-------------------------
xstep = datenum([0 0 1 0 0 0]);
fz = 14;
lw = 1.5;
% [r c]=find(isnan(tauW));          %È¥µôNAN
% tauC(r,:)= nan;
% tauCW_soulsby(r,:)= nan;
% tauCW_grant(r,:)= nan;
interval = datenum([0 0 0 0 30 0]);
blank1= find(diff(time1)>interval);
blank1 = [0;blank1;length(time1)];
%%

data2 = xlsread('d2_GT01_wave_parameters.xlsx');   
time2 = data2(:,1);            
tauW2= data2(:,14);
tauC2= data2(:,15);  %rs
tauCW_soulsby2= data2(:,22);
tauCW_grant2= data2(:,23);
blank2= find(diff(time2)>interval);
blank2 = [0;blank2;length(time2)];




%%
figure;
subplot(3,1,1);
gca1 = subplot(3,1,1);  
for ii = 1:length(blank1)-1
    index = blank1(ii)+1:blank1(ii+1);
    plot(time1(index),tauW1(index),'r.-'); hold on
end
hold on;
h1 = plot(time1(1),tauW1(1),'r.-');
for iii = 1:length(blank2)-1
    index = blank2(iii)+1:blank2(iii+1);
    plot(time2(index),tauW2(index),'b.-'); hold on
end
datetick('keepticks');
h2 = plot(time2(1),tauW2(1),'b.-');
%set(gca1,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca1,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\tauw (N/m^2)','Fontsize', 14);
legend([h1,h2],'CT01','GT01')
legend('boxoff');
set(gca1,'fontsize',14,'fontname','times new roman');
set(gca1,'linewidth',lw);
box off;
%%
subplot(3,1,2);
gca2=subplot(3,1,2);
for ii = 1:length(blank1)-1
    index = blank1(ii)+1:blank1(ii+1);
    plot(time1(index),tauC1(index),'r.-'); hold on
end
h1 = plot(time1(1),tauC1(1),'r.-');
hold on;
for iii = 1:length(blank2)-1
    index = blank2(iii)+1:blank2(iii+1);
    plot(time2(index),tauC2(index),'b.-'); hold on
end
h2 = plot(time2(1),tauC2(1),'b.-');
datetick('keepticks');
%set(gca2,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca2,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\tauc (N/m^2)','Fontsize', 14);
legend([h1,h2],'CT01','GT01')
legend('boxoff');
set(gca2,'fontsize',14,'fontname','times new roman');
set(gca2,'linewidth',lw);
box off;
%%
subplot(3,1,3);
gca3=subplot(3,1,3);
for ii = 1:length(blank1)-1
    index = blank1(ii)+1:blank1(ii+1);
    plot(time1(index),tauCW_soulsby1(index),'r.-'); hold on
end
hold on;
h1 = plot(time1(1),tauCW_soulsby1(1),'r.-');
for iii = 1:length(blank2)-1
    index = blank2(iii)+1:blank2(iii+1);
    plot(time2(index),tauCW_soulsby2(index),'b.-'); hold on
end
h2 = plot(time2(1),tauCW_soulsby2(1),'b.-');

datetick('keepticks');
%set(gca3,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca3,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\taucw (N/m^2)','Fontsize', 14);
legend([h1,h2],'CT01','GT01')
legend('boxoff');
set(gca3,'fontsize',14,'fontname','times new roman');
set(gca3,'linewidth',lw);
box off;