% 对adv数据进行burst平均并计算tauc
clear all; clc;
tic;
data1=load('H:\little_paper2\数据\N3_181225\adv4841\adv4841_despiking_majorminor_3beam.mat'); %水深
% data1.new_data1=data1.new_data2;
burst = data1.new_data1(:,1); time = data1.new_data1(:,2);
d = 0.696; %注意压力探头距离底部距离
Depth = data1.new_data1(:,3) + d; 
u_major = data1.new_data1(:,4); %主流向
v_minor = data1.new_data1(:,5); %次流向
w = data1.new_data1(:,6);
snr1 = data1.new_data1(:,7);snr2 = data1.new_data1(:,8);snr3 = data1.new_data1(:,9);
snr = (snr1+snr2+snr3)/3;
%
data2 = load('H:\little_paper2\数据\N3_181225\adv4841\adv4841_despiking_3beam.mat');
u = data2.new_data2(:,4); % 东向流速
v = data2.new_data2(:,5); % 北向流速
%u_major = v*cos(56.1157*pi/180)-u*cos((90-56.1157)*pi/180);
%v_minor = -v*sin(56.1157*pi/180)-u*sin((90-56.1157)*pi/180);
%% 构建数组存储变量
gap = find(diff(burst)~=0);
gap = [0;gap;length(u_major)-1];
for ii = 1:length(gap)-1
    index = gap(ii)+1:gap(ii+1);    % uu is u', uuww is u'w'
    DATA(ii)=struct('time',time(index),'burst',burst(index),'u_major',u_major(index),'v_minor',v_minor(index),'u',u(index),...
        'v',v(index),'w',w(index),'uu',u_major(index)-mean(u_major(index)),...
        'vv',v_minor(index)-mean(v_minor(index)),'ww',w(index)-mean(w(index)),'uu2',(u_major(index)-mean(u_major(index))).^2,...
        'vv2',(v_minor(index)-mean(v_minor(index))).^2,'ww2',(w(index)-mean(w(index))).^2,'depth',Depth(index),'snr',snr(index),...
        'uuww',(u_major(index)-mean(u_major(index))).*(w(index)-mean(w(index))));
end
%% 求burst平均的u_major,v_minor,u,v,w,uu,vv,ww,uu2,vv2,ww2,tke,mtke,cov
density = 1018;
for i = 1:length(gap)-1
    ave_burst(i,1)=mean(DATA(i).burst);
    ave_u_major(i,1)=mean(DATA(i).u_major);
    ave_v_minor(i,1)=mean(DATA(i).v_minor);
    ave_u(i,1)=mean(DATA(i).u);
    ave_v(i,1)=mean(DATA(i).v);
    ave_w(i,1)=mean(DATA(i).w);
    ave_uu(i,1) = mean(DATA(i).uu);  % u'
    ave_uu2(i,1) = mean(DATA(i).uu2);  % u'^2
    ave_vv(i,1) = mean(DATA(i).vv);
    ave_vv2(i,1) = mean(DATA(i).vv2);
    ave_ww(i,1) = mean(DATA(i).ww);
    ave_ww2(i,1) = mean(DATA(i).ww2);
    ave_uuww(i,1) = mean(DATA(i).uuww);
    ave_depth(i,1) = mean(DATA(i).depth);
    ave_snr(i,1) = mean(DATA(i).snr);
    ave_time(i,1) = mean(DATA(i).time);
    tauc_tke(i,1) = (ave_uu2(i,1)+ave_vv2(i,1)+ave_ww2(i,1))/2*0.2*density;  % TKE
    tauc_mtke(i,1) = density*0.9*ave_ww2(i,1); % M-TKE
    tauc_cov(i,1) = -density*ave_uuww(i,1); % RS
end
ave_U = sqrt(ave_u.^2+ave_v.^2); 
% ave_burstsave = [ave_burst ave_time ave_u_major ave_v_minor ave_w];
ave_burstsave = [ave_burst ave_time ave_u ave_v ave_w];
% ave_burstsave(:,2) = (ave_burstsave(:,2)-datenum([2018 7 13 19 0 0]))*24;
save('H:\little_paper2\数据\N3_181225\adv4841\adv4841_Ave_3beam.mat','-mat','ave_burstsave');
%% plot shear stress by current
figure; fz = 15;
plot(ave_time,abs(tauc_cov),'r'); hold on;
plot(ave_time,tauc_tke,'g');
plot(ave_time,tauc_mtke,'b');
datetick;
legend('COV','TKE','M-TKE'); legend('boxoff')
xlabel('Time(Month/Day)','fontsize',fz,'fontname','times new roman');
ylabel('Shear stess(N/m^2)','fontsize',fz,'fontname','times new roman');
title('Shear stress due to current','fontsize',fz,'fontname','times new roman');
set(gca,'fontsize',fz,'fontname','times new roman');
%%
str_time = datestr(ave_time); yr = year(str_time); mh = month(str_time); dy = day(str_time);
hr = hour(str_time); mt = minute(str_time); sd = second(str_time); 
parameters=[ave_burst yr mh dy hr mt sd ave_depth ave_U ave_u ave_v ave_w ave_snr];
xlswrite('burst_averaged_d1.xlsx',{'burst','year','month','day','hour','minute','second','Depth(m)','U(m/s)','u_major(m/s)',...
    'v_minor(m/s)','w(m/s)','ave_SNR'},'sheet1','a1');
xlswrite('burst_averaged_d1.xlsx',parameters,'sheet1','a2');
%% 潮流椭圆
figure(1)
a = max(ave_u_major); b = max(ave_v_minor); x0 = 0; y0 = 0; color = 'k'; 
u_tmp = [u v];
[theta,d,zr] = paxis(u_tmp);
phi = theta;
ellipse(a,b,phi,x0,y0,color)
axis equal
xlabel('East velocity')
ylabel('North velocity');
title('Tidal ellipse')
%%
a = 0.5999; b = 0.0099;


u =  lsqcurvefit(@(a,phi,t)a*cos(w*t-phi))
%%
%对东向和北向流速进行cosine函数拟合，计算潮流椭圆的长半轴和短半轴
p_u = fittype('a_u*cos(2*pi/12.46*24*ave_time-phi_u)','independent','ave_time');
f_u = fit(ave_time,ave_u,p_u);
plot(f_u,ave_time,ave_u);
a_u = f_u.a_u;
phi_u = f_u.phi_u*180/pi; %degree

figure(2)
p_v = fittype('a_v*cos(2*pi/12.46*24*ave_time-phi_v)','independent','ave_time');
f_v = fit(ave_time,ave_v,p_v);
plot(f_v,ave_time,ave_v);
a_v = f_v.a_v;
phi_v = f_v.phi_v*180/pi; %degree
[SEMA,ECC, INC,PHA]=ap2ep(a_u,phi_u,a_v,phi_v);
%%
figure(3)
ellipse(SEMA,SEMA*ECC,phi,x0,y0,color);
%%
% 将流速burst平均之后再进行cosine拟合
% http://blog.csdn.net/abcjennifer/article/details/7684836
% phi = theta-pi/2;
% plot(u_despiking,v_despiking,'.'); hold on


ellipse(a,b,phi,x0,y0,color)
axis equal
xlabel('East velocity')
ylabel('North velocity');
title('Tidal ellipse')

%% current direction   这种方法只针对 u,v 表示东向和北向流速
% dir = atan(ave_v./ave_u)*180/pi;
% for i = 1:length(ave_u)
%     if ave_u(i)>0&&ave_v(i)>0
%         dir(i) = 90 - dir(i) ;
%     else if ave_u(i)>0&&ave_v(i)<0
%             dir(i) = 90 - dir(i);
%         else if ave_u(i)<0&&ave_v(i)<0
%                 dir(i) = 270 - dir(i);
%             else if ave_u(i)<0&&ave_v(i)>0
%                     dir(i) = 270 - dir(i);
%                 end
%             end
%         end
%     end
% end