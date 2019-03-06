%% Procedure 计算波浪切应力的程序
function [touWave,Uwave]=touWaveProcedure(Hs,Tz,height,d50,density,viscosity)
%Hs 有效波高
%Tz 有效波高对应的波周期
%height 距底高度
%d50 中值粒径
%先计算光滑底床情况时摩擦系数fws，再根据计算的雷诺数判断是层流还是紊流。
%计算粗糙底床的摩擦系数fwr，如果fwr>fws，则选取粗糙的摩擦系数，如果粗糙系数还小于光滑时的摩擦系数，则选取光滑时的摩擦系数。

density=1027;%kg·m-3
viscosity=1.36*10^-6;%m2·s-1

H=Hs/sqrt(2);%m
Tp=Tz;%sT   %  CT02改过了  Tp=1.28*Tz
k=wavek(1/Tp,height);%无量纲
Uwave=pi*H/(Tp*sinh(k*height));%m·s-1

A=Uwave*Tp/(2*pi);%m
ReynoldsNumber=Uwave*A/viscosity;%无量纲
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
 fwr=0.237*r^-0.52;%无量纲
 if fwr>fws
     touWave=1/2*density*fwr*Uwave^2;
 else 
     touWave=1/2*density*fws*Uwave^2;
 end

