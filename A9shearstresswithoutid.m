clear all;
clc;
data = xlsread('adv_wave_parameters.xlsx');   
time = data(:,1);            
tauW= data(:,14);
tauC_rs= data(:,15);  %rs
tauC_tke= data(:,16);
tauC_mtke= data(:,17);
tauC_id= data(:,18);
tauCW_soulsby_rs= data(:,22);
tauCW_grant_rs= data(:,23);
tauCW_soulsby_tke= data(:,25);
tauCW_grant_tke= data(:,26);
tauCW_soulsby_mtke= data(:,28);
tauCW_grant_mtke= data(:,29);
tauCW_soulsby_id= data(:,31);
tauCW_grant_id= data(:,32);


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
blank= find(diff(time)>interval);
blank = [0;blank;length(time)];
figure;
subplot(3,1,1);
gca1 = subplot(3,1,1);  %
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauW(index),'b.-'); hold on
end
datetick('keepticks');
%set(gca1,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca1,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\tauw (N/m^2)','Fontsize', 14);
legend('¦Ów');
legend('boxoff');
set(gca1,'fontsize',14,'fontname','times new roman');
set(gca1,'linewidth',lw);
box off;
%%
subplot(3,1,2);
gca2=subplot(3,1,2);
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauC_rs(index),'b.-'); hold on
end
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauC_tke(index),'r.-'); hold on
end
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauC_mtke(index),'k.-'); hold on
end
% for ii = 1:length(blank)-1
%     index = blank(ii)+1:blank(ii+1);
%     plot(time(index),tauC_id(index),'g.-'); hold on
% end
A1=plot(time(1),tauC_rs(1),'b.-');
A2= plot(time(1),tauC_tke(1),'r.-');
A3=plot(time(1),tauC_mtke(1),'k.-');
datetick('keepticks');
%set(gca2,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca2,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\tauc (N/m^2)','Fontsize', 14);
legend([A1 A2 A3],'rs','tke','mtke');
legend('boxoff');
set(gca2,'fontsize',14,'fontname','times new roman');
set(gca2,'linewidth',lw);
box off;
%%
subplot(3,1,3);
gca3=subplot(3,1,3);
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauCW_soulsby_rs(index),'b.-'); hold on
end
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauCW_soulsby_tke(index),'r.-'); hold on
end
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauCW_soulsby_mtke(index),'k.-'); hold on
end
% for ii = 1:length(blank)-1
%     index = blank(ii)+1:blank(ii+1);
%     plot(time(index),tauCW_soulsby_id(index),'g.-'); hold on
% end
A1=plot(time(1),tauCW_soulsby_rs(1),'b.-');
A2= plot(time(1),tauCW_soulsby_tke(1),'r.-');
A3=plot(time(1),tauCW_soulsby_mtke(1),'k.-');
%A4= plot(time(1),tauCW_soulsby_id(1),'g.-')
datetick('keepticks');
%set(gca3,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca3,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\taucw (N/m^2)','Fontsize', 14);
legend([A1 A2 A3],'rs','tke','mtke');
legend('boxoff');
set(gca3,'fontsize',14,'fontname','times new roman');
set(gca3,'linewidth',lw);
box off;
%%
figure;
gca4 = subplot('position',[0.1 0.1 0.8 0.35]);

    for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauCW_grant_rs(index),'b.-'); hold on
end
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauCW_grant_tke(index),'r.-'); hold on
end
for ii = 1:length(blank)-1
    index = blank(ii)+1:blank(ii+1);
    plot(time(index),tauCW_grant_mtke(index),'k.-'); hold on
end
% for ii = 1:length(blank)-1
%     index = blank(ii)+1:blank(ii+1);
%     plot(time(index),tauCW_grant_id(index),'g.-'); hold on
% end
A1=plot(time(1),tauCW_grant_rs(1),'b.-');
A2= plot(time(1),tauCW_grant_tke(1),'r.-');
A3=plot(time(1),tauCW_grant_mtke(1),'k.-');
%A4= plot(time(1),tauCW_grant_id(1),'g.-')
datetick('keepticks');
%set(gca3,'xlim',[xmin xmax],'xtick',[xmin:xstep:xmax]);
%set(gca1,'ylim',[0 8000],'ytick',[0:4000:8000],'fontsize',fz,'fontname','times new roman');
set(gca4,'fontsize',fz,'fontname','times new roman');
datetick('keepticks');
xlabel('Time','Fontsize', 14,'fontname','times new roman');
ylabel('\taucw (N/m^2)','Fontsize', 14);
legend([A1 A2 A3],'rs','tke','mtke');
legend('boxoff');
set(gca4,'fontsize',14,'fontname','times new roman');
set(gca4,'linewidth',lw);
box off;