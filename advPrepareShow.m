% Despiking Acoustic Doppler Velocimeter Data  --Derek G.Goring & Viadimir I.Nikora
% function [u]=advPrepare(u);
function U=advPrepareShow(u)
c1=1.5;c2=1.35;
spike = [];
for i = 1:length(u)
    if i ~= 1 && i ~= length(u);
        deltaU1(i,1) = (u(i+1,1)-u(i-1,1))/2;
    else
        deltaU1(i,1)=(u(i,1)-mean(u,1))/2;
    end
end
for i=1:length(deltaU1)
    if i~=1&&i~=length(deltaU1);
        deltaU2(i,1)=(deltaU1(i+1,1)-deltaU1(i-1,1))/2;
    else
        deltaU2(i,1)=(deltaU1(i,1)-mean(deltaU1,1))/2;
    end
end
thetaU=median(abs(u-median(u)));
thetaDu=median(abs(deltaU1-median(deltaU1)));
thetaD2u=median(abs(deltaU2-median(deltaU2)));
URange=[-c1*thetaU+mean(u) c1*thetaU+mean(u)];
deltaURange=[-c1*thetaDu+mean(deltaU1,1) c1*thetaDu+mean(deltaU1,1)];
deltaU2Range=[-c1*thetaD2u+mean(deltaU2,1) c1*thetaD2u+mean(deltaU2,1)];
um=c2*thetaU*sqrt(2*log(length(u)));
deltaUm=c2*thetaDu*sqrt(2*log(length(deltaU1)));
deltaU2m=c2*thetaD2u*sqrt(2*log(length(deltaU2)));
thetaD2UU=atan(sum(u.*deltaU2)/sum(u.^2));

%%
duua=um;%u--x
duub=deltaUm;
d2udua=deltaUm;%du--x
d2udub=deltaU2m;
A=[(cos(thetaD2UU))^2 (sin(thetaD2UU))^2;(sin(thetaD2UU))^2 (cos(thetaD2UU))^2];
B=[um^2;deltaU2m^2];
X=A\B;
d2uua=sqrt(X(1));%u--x
d2uub=sqrt(X(2));
theta=0:pi/30:2*pi;

%% detect du_u
X1=duua.*cos(theta)+median(u);
Y1=duub.*sin(theta)+median(deltaU1);
figure;
plot(u,deltaU1,'*');
hold on
plot(X1,Y1,'r');
spikeTemp=find( u>max(X1)|u<min(X1));
spike=[spike;spikeTemp];
plot(u(spikeTemp),deltaU1(spikeTemp),'ro');
alpha=acos((u-median(u))./duua);
ugood=duub.*sin(alpha);
spikeTemp=find(abs(deltaU1)>ugood);
spikeTemp(find((u(spikeTemp)>=URange(1)&u(spikeTemp)<=URange(2))&(deltaU1(spikeTemp)>=deltaURange(1)&deltaU1(spikeTemp)<=deltaURange(2))))=[];
spike=[spike;spikeTemp];
plot(u(spikeTemp),deltaU1(spikeTemp),'ro');
%% detect d2u_du
X2=d2udua.*cos(theta)+median(deltaU1);
Y2=d2udub.*sin(theta)+median(deltaU2);
figure;
plot(deltaU1,deltaU2,'*');
hold on
plot(X2,Y2,'r');
spikeTemp=find( deltaU1>max(X2)|deltaU1<min(X2));
spike=[spike;spikeTemp];
plot(deltaU1(spikeTemp),deltaU2(spikeTemp),'ro');
alpha=acos((deltaU1-median(deltaU1))./d2udua);
ugood=d2udub.*sin(alpha);
spikeTemp=find(abs(deltaU2)>ugood);
spikeTemp(find((deltaU1(spikeTemp)>=deltaURange(1)&deltaU1(spikeTemp)<=deltaURange(2))&(deltaU2(spikeTemp)>=deltaU2Range(1)&deltaU2(spikeTemp)<=deltaU2Range(2))))=[];
spike=[spike;spikeTemp];
plot(deltaU1(spikeTemp),deltaU2(spikeTemp),'ro');
%% detect d2u_u,先平移，后旋转
X3=d2uua.*cos(theta).*cos(-thetaD2UU)+d2uub.*sin(theta).*sin(-thetaD2UU)+median(u);
Y3=-d2uua.*cos(theta).*sin(-thetaD2UU)+d2uub.*sin(theta).*cos(-thetaD2UU)+median(deltaU2);
figure;
plot(u,deltaU2,'*');
hold on
plot(X3,Y3,'r');
uR=(u-median(u)).*cos(thetaD2UU)+(deltaU2-median(deltaU2)).*sin(thetaD2UU);
deltaU2R=-(u-median(u)).*sin(thetaD2UU)+(deltaU2-median(deltaU2)).*cos(thetaD2UU);
spikeTemp=find( uR>max(d2uua.*cos(theta))|uR<min(d2uua.*cos(theta)));
spike=[spike;spikeTemp];
plot(u(spikeTemp),deltaU2(spikeTemp),'ro');
alpha=acos(uR./d2uua);
ugood=d2uub.*sin(alpha);
spikeTemp=find(abs(deltaU2R)>ugood);
spikeTemp(find((u(spikeTemp)>=URange(1)&u(spikeTemp)<=URange(2))&(deltaU2(spikeTemp)>=deltaU2Range(1)&deltaU2(spikeTemp)<=deltaU2Range(2))))=[];
spike=[spike;spikeTemp];
plot(u(spikeTemp),deltaU2(spikeTemp),'ro');
spike=unique(spike);
figure;
plot3(u,deltaU1,deltaU2,'.');
%% replace

u(spike)=NaN;%%make data NaN
u=replaceInf(u);%%replace it
if ~isempty(spike)
for i=1:length(spike)
    if spike(i)>12&&spike(i)+12<=length(u)
        xx=(1:25)';
        yy=u(spike(i)-12:spike(i)+12,1);
        [p,s]=polyfit(xx,yy,3);
        u(spike(i),1)=p(1)*13^3+p(2)*13^2+p(3)*13+p(4);
    elseif spike(i)<=12
        xx=(1:spike(i)+12)';
        yy=u(1:spike(i)+12,1);        
        [p,s]=polyfit(xx,yy,3);
        u(spike(i),1)=p(1)*spike(i)^3+p(2)*spike(i)^2+p(3)*spike(i)+p(4);
    elseif spike(i)+12>length(u)
        xx=(spike(i)-12:length(u))';
        yy=u(spike(i)-12:length(u),1);
        u(spike(i),1)=p(1)*13^3+p(2)*13^2+p(3)*13+p(4);
    end
end
end
U=u;