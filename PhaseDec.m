function [newU]=PhaseDec(RawU,Fs)
%phase method wave-turbulence decomposition
%RawU :raw data
%N :number of data points
%Fs :sampling frequence
%low frequence cut off=1/15;
%high frequence cut off=1;
hf=1/3;
lf=1/15;
N=length(RawU);
if mod(N,2)==1
    N=N+1;
end
[Spectral,F]=pwelch(RawU-mean(RawU),[],[],N,Fs,'twoside');%Welch ��������
F(1)=10^-16;
FouU=fft(RawU,N);%����Ҷ�任
SpeTemp=Spectral(1:N/2);
FTemp=F(1:N/2);
SpeTemp(F>=lf&F<=hf)=NaN;%����
FTemp(F>=lf&F<=hf)=NaN;
x=FTemp(~isnan(FTemp));
x(1)=10^-16;
y=SpeTemp(~isnan(SpeTemp));
yy=interp1(log10(x(1:end)),log10(y(1:end)),log10(F));%��ֵ
yy=10.^yy;
F(1)=0;
%Y=sqrt(yy(isnan(FTemp))*pi*N*Fs);%��ֵ
Y=yy(isnan(FTemp))*N*Fs;
FourTempL=FouU(isnan(FTemp));%��Ҫ�ı�ĸ���Ҷ��,���
FourTempR=FouU(F>=Fs-hf&F<=Fs-lf);%��Ҫ�ı�ĸ���Ҷ��Ҳ�
%phaseL=atan(imag(FourTempL)./real(FourTempL));
%phaseR=atan(imag(FourTempR)./real(FourTempR));
%NewFourL=sqrt(Y).*exp(1i*phaseL);
%NewFourR=sqrt(Y).*exp(1i*phaseR);
NewFourL=sqrt(Y)./abs(FourTempL).*FourTempL;
NewFourR=flipud(sqrt(Y))./abs(FourTempR).*FourTempR;
FouU(F>=Fs-hf&F<=Fs-lf)=NewFourR;
FouU(isnan(FTemp))=NewFourL;%��������Ҷ��
newU=real(ifft(FouU,N));

