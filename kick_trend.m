% We all know that spectral estimation request that the sigal is at least weakly stationary. So, detrending 
% is just trying to remove the slow nonstationary trends from the data that need processing.
% This is a detrending method using smoothness priors aproach(SPA).
% Ref. An advanced detrending method with application to HRV analysis. IEEE trans. on Biomedical Engineering

% An example for detrending. the sinusoidal signal that have lower frequency is used as the trend.
% n = 0:1/40:(1024-1)/40;
% data = sin(2*pi*0.1*n)+sin(2*pi*5*n);
% data = data(:);
% N = length(data);
% lambda = 100;
% I = speye(N);  % 单位稀疏矩阵
% D2 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
% trend = inv(I+lambda^2*D2'*D2)*data;
% detrenddata = data-trend;
% subplot(211);
% plot(n,data,'r',n,trend,'g');
% title('the orginal data and trend');
% legend('the orginal data','the trend');
% % xlim([0 5]);
% subplot(212);
% plot(n,detrenddata,'m')
% title('the data after detrenging');
% % xlim([0 5]);
%% This function is used for low frequency trend term elimination.去除低频信号
function [kt_data]=kick_trend(data,lambda)
% data is the original data with trend term
% lambda is a parameter that control the cutoff frequence, the
% corresponding cutoff frequencies is y=0.189*lambda^(-0.506)*sampling
% frequency(experiential).频率控制式子，低于该频率的信号将被剔除。

data = data(:);
N = length(data);
lambda = lambda;
I = speye(N);
D2 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);%N-2表示什么？如果你了解了，请告诉我。
trend = (I+lambda^2*(D2')*D2)\data;
kt_data = data-trend;
