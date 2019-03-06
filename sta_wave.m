function [Hm,Ha,H3,wT,HT,mT,H10,L]=sta_wave(p,dt,depth,hp,hf2,factor)
%ֱ��ͳ�Ƶõ����ֲ���
%Hm:��󲨸�
%Ha:ƽ������
%H3:1/3��󲨸�,����Ч����
%H10:1/10��󲨸ߣ���һ����Ч����
%wT����󲨸߶�Ӧ�Ĳ����ڣ�
%mT��ƽ����Ч������

mean_pressure=mean(p);
len_pressure=length(p);
index_0=0;
for index=1:len_pressure-1;
    if p(index)<=mean_pressure&&p(index+1)>mean_pressure
        index_0=index_0+1;
        wave_index(index_0,1)=index;%�ҵ�������λ��.ÿ���������֮���ʾһ������
    end
    
end
clear index_0
wave_time = diff(wave_index).*dt;
for index=1:length(wave_index)-1
    %wave_time(index,1)=(wave_index(index+1)-wave_index(index)+1)*dt;%ÿ�����εĲ�����
    wave_height(index,1)=max(p(wave_index(index):wave_index(index+1)))-min(p(wave_index(index):wave_index(index+1)));%ÿ�������ڵ����ֵ�����ֵ�Ĳ���ǲ���
end


    sigWaveTime=wave_time(wave_time>(1/hf2));%������Ƶ�ض�Ƶ�ʵĲ����������㿼�ǣ���Ϊ��Ч
    sigWaveHeight=wave_height(wave_time>(1/hf2));%����ʱ���жϵõ�����Ч�Ĳ���
if length(sigWaveTime)<=3;
   disp('there is no enough waves in this data section, please try another one');%������ĿС��3�����߱�ͳ������
   
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
        k=wavek(1/sigWaveTime(index),depth);%���㲨��
        Kz=(cosh(k*hp)/cosh(k*depth)).^1;
        L(index)=2*pi/k;
        if Kz<factor%����˥��ϵ������ֵ
            Kz=factor*2.5;%��Ϊ����һ��У��ϵ��
            sigWaveHeight(index)=sigWaveHeight(index)/Kz;
        else
            sigWaveHeight(index)=sigWaveHeight(index)/Kz;
        end
    end

[sigWaveHeight_temp,index]=sort(sigWaveHeight,'descend');%���������򣬴ӵ͵���
Ha=mean(sigWaveHeight);
HT=mean(sigWaveTime);
%L=mean(L);
H3=mean(sigWaveHeight_temp(1:fix(length(sigWaveHeight_temp)/3)));%1/3��Ч����
H10=mean(sigWaveHeight_temp(1:fix(length(sigWaveHeight_temp)/10)));%1/10��Ч����
mT=mean(sigWaveTime(index(1:fix(length(sigWaveHeight_temp)/3))));%1/3
[Hm,local]=max(sigWaveHeight);
wT=sigWaveTime(local);
L=L(local);
end
