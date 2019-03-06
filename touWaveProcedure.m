%% Procedure ���㲨����Ӧ���ĳ���
function [touWave,Uwave]=touWaveProcedure(Hs,Tz,height,d50,density,viscosity)
%Hs ��Ч����
%Tz ��Ч���߶�Ӧ�Ĳ�����
%height ��׸߶�
%d50 ��ֵ����
%�ȼ���⻬�״����ʱĦ��ϵ��fws���ٸ��ݼ������ŵ���ж��ǲ�������������
%����ֲڵ״���Ħ��ϵ��fwr�����fwr>fws����ѡȡ�ֲڵ�Ħ��ϵ��������ֲ�ϵ����С�ڹ⻬ʱ��Ħ��ϵ������ѡȡ�⻬ʱ��Ħ��ϵ����

density=1027;%kg��m-3
viscosity=1.36*10^-6;%m2��s-1

H=Hs/sqrt(2);%m
Tp=Tz;%sT   %  CT02�Ĺ���  Tp=1.28*Tz
k=wavek(1/Tp,height);%������
Uwave=pi*H/(Tp*sinh(k*height));%m��s-1

A=Uwave*Tp/(2*pi);%m
ReynoldsNumber=Uwave*A/viscosity;%������
if ReynoldsNumber<=5*10^5;
    B=2;N=0.5;
    fws=B*ReynoldsNumber^-N;
else ReynoldsNumber>5*10^5;
    B=0.0521;N=0.187;
    fws=B*ReynoldsNumber^-N;
end
%touWave=1/2*density*fws*Uwave^2;
 
 ks=2.5*d50;%mm
 %ks=0.041;%m
 r=A/ks;%
 fwr=0.237*r^-0.52;%������
 if fwr>fws
     touWave=1/2*density*fwr*Uwave^2;
 else 
     touWave=1/2*density*fws*Uwave^2;
 end

