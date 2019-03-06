%-----------------ADV数据处理第二步:毛刺点的去除与替换------------------------
tic;
clc;clear all;close all;
%
DATA = load('CT01_02_qualified.mat');
data1= DATA.data(1:200000,:);
data2= DATA.data(200001:239730,:);
data6= DATA.data(240000:end,:);
data= DATA.data(239730:239999,:);
burst = data(:,1); 
u = data(:,4); v = data(:,5); w = data(:,6);
tol=10^-6;
gordan1=find(abs(u)>abs((10*mean(u))));
u(gordan1,:)=(u(gordan1+1,:)+u(gordan1-1,:))/2;
gordan2=find(abs(v)>abs((10*mean(v))));
v(gordan2,:)=(v(gordan2+1,:)+v(gordan2-1,:))/2;
gordan3=find(abs(w)>abs((10*mean(w))));
w(gordan3,:)=(w(gordan3+1,:)+w(gordan3-1,:))/2;
u_despiking=u;
v_despiking=v;
w_despiking=w;
% clear u1 v1 w1 tol
%% 去毛刺前后概率密度曲线
figure
subplot(2,2,1);
[f_u1,x_u1] = ksdensity(u);
plot(x_u1,f_u1,'b*'); hold on
[f_u2, x_u2] = ksdensity(u_despiking);
plot(x_u2,f_u2,'ro')
xlabel('East velocity (m/s)');
ylabel('PDF');
a_u = 0; sigma_u = sqrt(var(u)); % 高斯分布
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
a_v = 0; sigma_v = sqrt(var(v)); % 高斯分布
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
a_w = 0; sigma_w = sqrt(var(w)); % 高斯分布
x_w = min(w):0.0001:max(w);
y_w = (1/((sqrt(2*pi))*sigma_w))*exp(-((x_w-a_w).^2)/(2*sigma_w.^2));
plot(x_w,y_w,'k','LineWidth',1.5);
legend('w','w-despiking','Normal distribution');
clear f_w1 x_w1 f_w2 x_w2 a_w sigma_w x_w y_w

%% 去毛刺前后结果对比
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
%% 将流速分解到主流向和次流向上
%data = load('59484112_edu_despiking.mat');
% u=mu;v=mv;
%u_despiking = data.new_data2(:,4); v_despiking = data.new_data2(:,5); 

u_tmp = [u_despiking; v_despiking]';
[theta,d,zr] = paxis(u_tmp);
u_major = zr(:,1); %主流向
v_minor = zr(:,2); %次流向
% save
new_data1 = DATA.data;
new_data1(:,4) = u_major;
new_data1(:,5) = v_minor;
new_data1(:,6) = w_despiking;
save('CT01_02_despiking_major_minor.mat','new_data1');
new_data2 = DATA.data;
new_data2(:,4) = u_despiking;
new_data2(:,5) = v_despiking;
new_data2(:,6) = w_despiking;
save('CT01_02_qualified_despiking.mat','new_data2');



