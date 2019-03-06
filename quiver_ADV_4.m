tic;
clear all;close all;
clc;
DATA = load('CT02_despiking_aveUVW.mat');
Tm = DATA.aveUVW(:,1);
Um = DATA.aveUVW(:,2);
Vm = DATA.aveUVW(:,3);
%%
% windspeed = num_1(:,2); direction = num_1(:,3);
% hour1 = num_1(:,6); e_wind = num_1(:,4); n_wind = num_1(:,5);
y = ones(1,length(Tm))-1;
%
fz = 15; lw = 1.5;
figure;
subplot(3,1,1);
N=3;
Tm1 = Tm(1:N:end); y1 = y(1:N:end);Um1 = Um(1:N:end); Vm1 = Vm(1:N:end);
quiver(Tm1,y1,Um1,Vm1,'color','k','maxheadsize',0.01,'linewidth',1,'AutoScale','off');
hold on
plot(Tm1,y1,'k-','linewidth',lw);
datetick('keepticks');
set(gca,'linewidth',lw);
%set(gca1,'xlim',[12 210],'xtick',[12:12:210],'xticklabel',[],'linewidth',lw);
set(gca,'ylim',[-0.04 0.04],'ytick',[-0.04:0.02:0.04],'fontsize',fz,'fontname','times new roman');
ylabel('Flow velocity (m/s)','fontsize',fz,'fontname','times new roman');
box off;
legend('CT01-5195');
h1 = legend ('CT01-5195')
set (h1,'Fontsize', 10);%
