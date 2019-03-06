clc;clear all;close all;
%
DATA = load('advClean.mat');
u = DATA.data(:,4); v = DATA.data(:,5); w = DATA.data(:,6);
tol=10^-6;
%%
gap = find(diff(burst)~=0); gap = [0;gap;length(burst)];
% save spike data propotion
sper1(1:length(gap)-1) = 0.0; sper2(1:length(gap)-1) = 0.0; sper3(1:length(gap)-1) = 0.0;
for index = 1:length(gap)-1
    %index = 1:length(gap)-1
    uTemp = u(gap(index)+1:gap(index+1));  %u
    u1 = uTemp;
    [uTemp,spike] = advPrepare(uTemp);
    sper1(index) = length(spike)/length(u1)*100;
    while abs(u1-uTemp) >= tol;
        u1 = uTemp;
        uTemp = advPrepare(uTemp);
    end
    u_despiking(gap(index)+1:gap(index+1)) = uTemp;
    
    vTemp = v(gap(index)+1:gap(index+1));  %v
    v1 = vTemp;
    [vTemp,spike] = advPrepare(vTemp);
    sper2(index) = length(spike)/length(v1)*100;
    while abs(v1-vTemp) >= tol;
        v1 = vTemp;
        vTemp = advPrepare(vTemp);
    end
    v_despiking(gap(index)+1:gap(index+1)) = vTemp;
    
    wTemp = w(gap(index)+1:gap(index+1)); %w
    w1 = wTemp;
    [wTemp,spike] = advPrepare(wTemp);
    sper3(index) = length(spike)/length(w1)*100;
    while abs(w1-wTemp) >= tol;
        w1 = wTemp;
        wTemp = advPrepare(wTemp);
    end
    w_despiking(gap(index)+1:gap(index+1)) = wTemp;
end
% clear u1 v1 w1 tol
%% plot density distribution
figure
subplot(2,2,1);
[f_u1,x_u1] = ksdensity(u);
plot(x_u1,f_u1,'b*'); hold on
[f_u2, x_u2] = ksdensity(u_despiking);
plot(x_u2,f_u2,'ro')
xlabel('East velocity (m/s)');
ylabel('PDF');
a_u = 0; sigma_u = sqrt(var(u)); 
x_u = min(u):0.0001:max(u);
y_u = (1/((sqrt(2*pi))*sigma_u))*exp(-((x_u-a_u).^2)/(2*sigma_u.^2));
plot(x_u,y_u,'k','LineWidth',1.5);
legend('u','u-despiking','Normal distribution');
clear f_u1 x_u1 f_u2 x_u2 a_u sigma_u x_u y_u
%
subplot(2,2,2);
[f_v1,x_v1] = ksdensity(v);
plot(x_v1,f_v1,'b*'); hold on
[f_v2, x_v2] = ksdensity(v_despiking);
plot(x_v2,f_v2,'ro')
xlabel('North velocity (m/s)');
ylabel('PDF');
a_v = 0; sigma_v = sqrt(var(v));
x_v = min(v):0.0001:max(v);
y_v = (1/((sqrt(2*pi))*sigma_v))*exp(-((x_v-a_v).^2)/(2*sigma_v.^2));
plot(x_v,y_v,'k','LineWidth',1.5);
legend('v','v-despiking','Normal distribution');
clear f_v1 x_v1 f_v2 x_v2 a_v sigma_v x_v y_v

%
subplot(2,2,3);
[f_w1,x_w1] = ksdensity(w);
plot(x_w1,f_w1,'b*'); hold on
[f_w2, x_w2] = ksdensity(w_despiking);
plot(x_w2,f_w2,'ro')
xlabel('Vertical velocity (m/s)');
ylabel('PDF');
a_w = 0; sigma_w = sqrt(var(w));
x_w = min(w):0.0001:max(w);
y_w = (1/((sqrt(2*pi))*sigma_w))*exp(-((x_w-a_w).^2)/(2*sigma_w.^2));
plot(x_w,y_w,'k','LineWidth',1.5);
legend('w','w-despiking','Normal distribution');
clear f_w1 x_w1 f_w2 x_w2 a_w sigma_w x_w y_w

%% compare raw data and data after despiking
figure
subplot(3,1,1)
plot(u,'b.-'); hold on; plot(u_despiking,'r.-');
legend('original u','despiking u');
subplot(3,1,2)
plot(v,'b.-'); hold on; plot(v_despiking,'r.-');
legend('original v','despiking v');
subplot(3,1,3)
plot(w,'b.-'); hold on; plot(w_despiking,'r.-');
legend('original w','despiking w');
%
%% transform uEast, vNorth to uMajor, vMinor
%data = load('59484112_edu_despiking.mat');
% u=mu;v=mv;
%u_despiking = data.new_data2(:,4); v_despiking = data.new_data2(:,5); 

u_tmp = [u_despiking; v_despiking]';
[theta,d,zr] = paxis(u_tmp);
u_major = zr(:,1); % uMajor
v_minor = zr(:,2); % vMinor
% save
new_data1 = DATA.data;
new_data1(:,4) = u_major;
new_data1(:,5) = v_minor;
new_data1(:,6) = w_despiking;
save('advClean_despiking_majorminor.mat','new_data1');
new_data2 = DATA.data;
new_data2(:,4) = u_despiking;
new_data2(:,5) = v_despiking;
new_data2(:,6) = w_despiking;
save('advClean_despiking.mat','new_data2');

%% plot tidal ellipse
figure
a = d(1)*mean(u_despiking); b = d(2)*mean(v_despiking); x0 = 0; y0 = 0; color = 'k'; 
phi = theta;
% phi = theta-pi/2;
% plot(u_despiking,v_despiking,'.'); hold on
ellipse(a,b,phi,x0,y0,color)
axis equal
xlabel('East velocity')
ylabel('North velocity');
title('Tidal ellipse')
toc;


