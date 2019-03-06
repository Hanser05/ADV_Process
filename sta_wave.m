function [Hm,Ha,H3,wT,HT,mT,H10,L]=sta_wave(p,dt,depth,hp,hf2,factor)
%直接统计得到各种波高
%Hm:最大波高
%Ha:平均波高
%H3:1/3最大波高,即有效波高
%H10:1/10最大波高，另一种有效波高
%wT：最大波高对应的波周期，
%mT：平均有效波周期

mean_pressure=mean(p);
len_pressure=length(p);
index_0=0;
for index=1:len_pressure-1;
    if p(index)<=mean_pressure&&p(index+1)>mean_pressure
        index_0=index_0+1;
        wave_index(index_0,1)=index;%找到跨零点的位置.每两个跨零点之间表示一个波形
    end
    
end
clear index_0
wave_time = diff(wave_index).*dt;
for index=1:length(wave_index)-1
    %wave_time(index,1)=(wave_index(index+1)-wave_index(index)+1)*dt;%每个波形的波周期
    wave_height(index,1)=max(p(wave_index(index):wave_index(index+1)))-min(p(wave_index(index):wave_index(index+1)));%每个波形内的最高值和最低值的差，就是波高
end


    sigWaveTime=wave_time(wave_time>(1/hf2));%超过高频截断频率的波，不作计算考虑，认为无效
    sigWaveHeight=wave_height(wave_time>(1/hf2));%根据时间判断得到的有效的波高
if length(sigWaveTime)<=3;
   disp('there is no enough waves in this data section, please try another one');%波浪数目小于3个不具备统计意义
   
   Hm=NaN;
   wT=NaN;
   Ha=NaN;
   mT=NaN;
   H3=NaN;
   H10=NaN;
   HT=NaN;
   L=NaN;
else 
    for index=1:length(sigWaveTime)
        k=wavek(1/sigWaveTime(index),depth);%计算波数
        Kz=(cosh(k*hp)/cosh(k*depth)).^1;
        L(index)=2*pi/k;
        if Kz<factor%波浪衰减系数的阈值
            Kz=factor*2.5;%人为赋予一个校正系数
            sigWaveHeight(index)=sigWaveHeight(index)/Kz;
        else
            sigWaveHeight(index)=sigWaveHeight(index)/Kz;
        end
    end

[sigWaveHeight_temp,index]=sort(sigWaveHeight,'descend');%将波高排序，从低到高
Ha=mean(sigWaveHeight);
HT=mean(sigWaveTime);
%L=mean(L);
H3=mean(sigWaveHeight_temp(1:fix(length(sigWaveHeight_temp)/3)));%1/3有效波高
H10=mean(sigWaveHeight_temp(1:fix(length(sigWaveHeight_temp)/10)));%1/10有效波高
mT=mean(sigWaveTime(index(1:fix(length(sigWaveHeight_temp)/3))));%1/3
[Hm,local]=max(sigWaveHeight);
wT=sigWaveTime(local);
L=L(local);
end
