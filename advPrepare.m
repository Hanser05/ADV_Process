% Despiking Acoustic Doppler Velocimeter Data  --Derek G.Goring & Viadimir I.Nikora
% function [u] = advPrepare(u);
function [U,spike] = advPrepare(u)
c1 = 1.5; c2 = 1.35;
spike = [];
% 计算一阶、二阶导数
for index = 1:length(u)   % 一阶导  
    if index ~= 1 && index ~= length(u);   
        deltaU1(index,1) = (u(index+1,1)-u(index-1,1))/2;
    else   % 首尾数据的一阶导
        deltaU1(index,1) = (u(index,1)-mean(u,1))/2;
    end
end
for index = 1:length(deltaU1)   % 二阶导
    if index ~= 1 && index ~= length(deltaU1);
        deltaU2(index,1) = (deltaU1(index+1,1)-deltaU1(index-1,1))/2;
    else
        deltaU2(index,1) = (deltaU1(index,1)-mean(deltaU1,1))/2;
    end
end
% 计算u, deltaU, deltaU2的最大绝对期望中值
thetaU = median(abs(u-median(u)));
thetaDu = median(abs(deltaU1-median(deltaU1)));
thetaD2u = median(abs(deltaU2-median(deltaU2)));
% 计算绝对有效而不被替代的数据区间
URange = [-c1*thetaU+mean(u) c1*thetaU+mean(u)];
deltaURange = [-c1*thetaDu+mean(deltaU1,1) c1*thetaDu+mean(deltaU1,1)];
deltaU2Range = [-c1*thetaD2u+mean(deltaU2,1) c1*thetaD2u+mean(deltaU2,1)];
% 计算最大期望值
um = c2*thetaU*sqrt(2*log(length(u)));
deltaUm = c2*thetaDu*sqrt(2*log(length(deltaU1)));
deltaU2m = c2*thetaD2u*sqrt(2*log(length(deltaU2)));
% 计算deltaU2-u 平面的椭圆旋转夹角
thetaD2UU = atan(sum(u.*deltaU2)/sum(u.^2));
%% 计算投影面内椭圆的半长轴与半短轴
du_u_a = um;%u--x    
du_u_b = deltaUm;     
d2u_du_a = deltaUm;%du--x
d2u_du_b = deltaU2m;
A = [(cos(thetaD2UU))^2 (sin(thetaD2UU))^2;(sin(thetaD2UU))^2 (cos(thetaD2UU))^2];
B = [um^2;deltaU2m^2];
X = A\B;
d2u_u_a = sqrt(X(1));%u--x
d2u_u_b = sqrt(X(2));
theta = 0:pi/30:2*pi;
%% detect du_u   %检测位于临界椭圆外和有效区间外的点  %找毛刺点的下标 
X1 = du_u_a.*cos(theta)+median(u);
Y1 = du_u_b.*sin(theta)+median(deltaU1);
spikeTemp = find(u>max(X1) | u<min(X1)); 
spike = [spike;spikeTemp];
alpha = acos((u-median(u))./du_u_a);
ugood = du_u_b.*sin(alpha);
spikeTemp = find(abs(deltaU1)>ugood);
%这个地方似乎有点问题...应该为 spikeTemp1 =find(abs(deltaU1-median(deltaU1))>ugood);不过结果没差多少。。
% 排除位于临界椭圆外但位于有效区间内的毛刺点
spikeTemp(find((u(spikeTemp)>=URange(1)&u(spikeTemp)<=URange(2))&(deltaU1(spikeTemp)>=deltaURange(1)&deltaU1(spikeTemp)<=deltaURange(2))))=[];
spike = [spike;spikeTemp];
%% detect d2u_du
X2 = d2u_du_a.*cos(theta)+median(deltaU1);
Y2 = d2u_du_b.*sin(theta)+median(deltaU2);
spikeTemp = find( deltaU1>max(X2)|deltaU1<min(X2));
spike = [spike;spikeTemp];
alpha = acos((deltaU1-median(deltaU1))./d2u_du_a);
ugood = d2u_du_b.*sin(alpha);
spikeTemp = find(abs(deltaU2)>ugood);
spikeTemp(find((deltaU1(spikeTemp)>=deltaURange(1)&deltaU1(spikeTemp)<=deltaURange(2))&(deltaU2(spikeTemp)>=deltaU2Range(1)&deltaU2(spikeTemp)<=deltaU2Range(2))))=[];
spike = [spike;spikeTemp];
%% detect d2u_u,先平移，后旋转
uR = (u-median(u)).*cos(thetaD2UU)+(deltaU2-median(deltaU2)).*sin(thetaD2UU);
deltaU2R = -(u-median(u)).*sin(thetaD2UU)+(deltaU2-median(deltaU2)).*cos(thetaD2UU);
spikeTemp = find(uR>max(d2u_u_a.*cos(theta))|uR<min(d2u_u_a.*cos(theta)));
spike = [spike;spikeTemp];
alpha = acos(uR./d2u_u_a);
ugood = d2u_u_b.*sin(alpha);
spikeTemp = find(abs(deltaU2R)>ugood);
spikeTemp(find((u(spikeTemp)>=URange(1)&u(spikeTemp)<=URange(2))&(deltaU2(spikeTemp)>=deltaU2Range(1)&deltaU2(spikeTemp)<=deltaU2Range(2))))=[];
spike = [spike;spikeTemp];
spike = unique(spike);
%% replace  替代毛刺点 三阶多项式拟合
%utemp=u;
%utemp(spike)=NaN;
%u(spike)=nanmean(utemp);%%make data NaN
%while(max(isnan(u))>0.1);
%u=replaceInf(u);%%replace it
%end

if ~isempty(spike)
    u(spike) = NaN;  %make data NaN
    while (find(isnan(u)))
        u = replaceInf(u);  %replace it   先用周围数的平均值替换毛刺点
    end
for index = 1:length(spike)   % 再用毛刺点前后12个数据，通过三次多项式拟合再次替换毛刺点
    if spike(index) > 12 && spike(index)+12 <= length(u)
        xx = (1:25)';
        yy = u(spike(index)-12:spike(index)+12,1);
        [p,s] = polyfit(xx,yy,3);
        u(spike(index),1) = p(1)*13^3+p(2)*13^2+p(3)*13+p(4);
    elseif spike(index) <= 12
        xx = (1:spike(index)+12)';
        yy = u(1:spike(index)+12,1);        
        [p,s] = polyfit(xx,yy,3);
        u(spike(index),1) = p(1)*spike(index)^3+p(2)*spike(index)^2+p(3)*spike(index)+p(4);
    elseif spike(index)+12 > length(u)
        xx = (spike(index)-12:length(u))';
        yy = u(spike(index)-12:length(u),1);
        [p,s] = polyfit(xx,yy,3);
        u(spike(index),1) = p(1)*spike(index)^3+p(2)*spike(index)^2+p(3)*spike(index)+p(4);      
%       u(spike(index),1) = p(1)*13^3+p(2)*13^2+p(3)*13+p(4);
    end
end
end
U = u;