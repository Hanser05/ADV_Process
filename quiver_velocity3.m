tic;
clear all;close all;
clc;
DATA = load('CT02_qualified_despiking_aveUVW.mat');
Tm = DATA.aveUVW(:,1);
Um = DATA.aveUVW(:,2);
Vm = DATA.aveUVW(:,3);
DATA2 = load('CT02_qualified_despiking_aveUVW.mat');
Tm2 = DATA2.aveUVW(:,1);
Um2 = DATA2.aveUVW(:,2);
Vm2 = DATA2.aveUVW(:,3);
%%
% windspeed = num_1(:,2); direction = num_1(:,3);
% hour1 = num_1(:,6); e_wind = num_1(:,4); n_wind = num_1(:,5);
y = ones(1,length(Tm))-1;
%
fz = 15; lw = 1.5;
figure
subplot(3,1,1)
gca1 = subplot(3,1,1);
% e_wind = e_wind(1:stride:end); n_wind = n_wind(1:stride:end);
N=3;
Tm1 = Tm(1:N:end); y1 = y(1:N:end);Um1 = Um(1:N:end); Vm1 = Vm(1:N:end);
quiver(Tm1,y1,Um1,Vm1,'color','k','maxheadsize',0.01,'linewidth',1,'AutoScale','off');
hold on
plot(Tm1,y1,'k-','linewidth',lw);
datetick('keepticks');
%set(gca,'xlim'[],'xtick','xticklabel','linewidth',lw);
set(gca1,'ylim',[-0.4 0.4],'ytick',[-0.4:0.2:0.4],'fontsize',fz,'fontname','times new roman');
ylabel('Flow velocity (m/s)','fontsize',fz,'fontname','times new roman');
box off;
legend('5194');
h1 = legend ('5194')
set (h1,'Fontsize', 8)
%
%%
y2 = ones(1,length(Tm2))-1;
subplot(3,1,2)
gca2 = subplot(3,1,2);
% e_wind = e_wind(1:stride:end); n_wind = n_wind(1:stride:end);
N=3;
Tm3 = Tm2(1:N:end); y3 = y2(1:N:end);Um3 = Um2(1:N:end); Vm3 = Vm2(1:N:end);
quiver(Tm3,y3,Um3,Vm3,'color','k','maxheadsize',0.01,'linewidth',1,'AutoScale','off');
hold on
plot(Tm3,y3,'k-','linewidth',lw);
datetick('keepticks');
% set(gca,'xlim','xtick','xticklabel','linewidth',lw);
set(gca2,'ylim',[-0.4 0.4],'ytick',[-0.4:0.2:0.4],'fontsize',fz,'fontname','times new roman');
ylabel('Flow velocity (m/s)','fontsize',fz,'fontname','times new roman');
box off
legend('4841')
h2 = legend ('4841')
set (h2, 'Fontsize', 8)
%%
subplot(3,1,3);
gca3 = subplot(3,1,3);
quiver(Tm1,y1,Um1,Vm1,'color','k','maxheadsize',0.01,'linewidth',1,'AutoScale','off');
hold on
datetick('keepticks');
%set(gca,'xlim'[],'xtick','xticklabel','linewidth',lw);
%set(gca3,'ylim',[-0.4 0.4],'ytick',[-0.4:0.2:0.4],'fontsize',fz,'fontname','times new roman');
ylabel('Flow velocity (m/s)','fontsize',fz,'fontname','times new roman');
hold on;
quiver(Tm3,y3,Um3,Vm3,'color','r','maxheadsize',0.01,'linewidth',1,'AutoScale','off');
hold on
%plot(Tm3,y3,'k-','linewidth',lw);
datetick('keepticks');
%set(gca3,'xlim','xtick','xticklabel','linewidth',lw);
set(gca3,'ylim',[-0.4 0.4],'ytick',[-0.4:0.2:0.4],'fontsize',fz,'fontname','times new roman');
ylabel('Flow velocity (m/s)','fontsize',fz,'fontname','times new roman');
box off
plot(Tm1,y1,'k-','linewidth',lw);
legend('5194','4841');
h = legend('5194','4841');
set(h,'Fontsize',8);