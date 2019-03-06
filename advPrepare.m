% Despiking Acoustic Doppler Velocimeter Data  --Derek G.Goring & Viadimir I.Nikora
% function [u] = advPrepare(u);
function [U,spike] = advPrepare(u)
c1 = 1.5; c2 = 1.35;
spike = [];
% ����һ�ס����׵���
for index = 1:length(u)   % һ�׵�  
    if index ~= 1 && index ~= length(u);   
        deltaU1(index,1) = (u(index+1,1)-u(index-1,1))/2;
    else   % ��β���ݵ�һ�׵�
        deltaU1(index,1) = (u(index,1)-mean(u,1))/2;
    end
end
for index = 1:length(deltaU1)   % ���׵�
    if index ~= 1 && index ~= length(deltaU1);
        deltaU2(index,1) = (deltaU1(index+1,1)-deltaU1(index-1,1))/2;
    else
        deltaU2(index,1) = (deltaU1(index,1)-mean(deltaU1,1))/2;
    end
end
% ����u, deltaU, deltaU2��������������ֵ
thetaU = median(abs(u-median(u)));
thetaDu = median(abs(deltaU1-median(deltaU1)));
thetaD2u = median(abs(deltaU2-median(deltaU2)));
% ���������Ч�������������������
URange = [-c1*thetaU+mean(u) c1*thetaU+mean(u)];
deltaURange = [-c1*thetaDu+mean(deltaU1,1) c1*thetaDu+mean(deltaU1,1)];
deltaU2Range = [-c1*thetaD2u+mean(deltaU2,1) c1*thetaD2u+mean(deltaU2,1)];
% �����������ֵ
um = c2*thetaU*sqrt(2*log(length(u)));
deltaUm = c2*thetaDu*sqrt(2*log(length(deltaU1)));
deltaU2m = c2*thetaD2u*sqrt(2*log(length(deltaU2)));
% ����deltaU2-u ƽ�����Բ��ת�н�
thetaD2UU = atan(sum(u.*deltaU2)/sum(u.^2));
%% ����ͶӰ������Բ�İ볤��������
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
%% detect du_u   %���λ���ٽ���Բ�����Ч������ĵ�  %��ë�̵���±� 
X1 = du_u_a.*cos(theta)+median(u);
Y1 = du_u_b.*sin(theta)+median(deltaU1);
spikeTemp = find(u>max(X1) | u<min(X1)); 
spike = [spike;spikeTemp];
alpha = acos((u-median(u))./du_u_a);
ugood = du_u_b.*sin(alpha);
spikeTemp = find(abs(deltaU1)>ugood);
%����ط��ƺ��е�����...Ӧ��Ϊ spikeTemp1 =find(abs(deltaU1-median(deltaU1))>ugood);�������û����١���
% �ų�λ���ٽ���Բ�⵫λ����Ч�����ڵ�ë�̵�
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
%% detect d2u_u,��ƽ�ƣ�����ת
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
%% replace  ���ë�̵� ���׶���ʽ���
%utemp=u;
%utemp(spike)=NaN;
%u(spike)=nanmean(utemp);%%make data NaN
%while(max(isnan(u))>0.1);
%u=replaceInf(u);%%replace it
%end

if ~isempty(spike)
    u(spike) = NaN;  %make data NaN
    while (find(isnan(u)))
        u = replaceInf(u);  %replace it   ������Χ����ƽ��ֵ�滻ë�̵�
    end
for index = 1:length(spike)   % ����ë�̵�ǰ��12�����ݣ�ͨ�����ζ���ʽ����ٴ��滻ë�̵�
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