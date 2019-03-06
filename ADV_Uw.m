function [Uw] = ADV_Uw(RawU,RawV,lf,hf,Fs)
% obtain Uw direcly from near bed velocity measurement by AVD
% INPUT PARAMETER
%         RawU: eastern velocity of each burst
%         RawV: northern velocity of each burst
%         lf: low frequency cutoff of wave
%         hf: high frequency cutoff of wave
%         Fs: smapling frequency
N = length(RawU);
[S_RawU,f] = pwelch(detrend(RawU),[],[],N,Fs,'oneside');
TurU = PhaseDec(RawU,Fs);
[S_turU,f] = pwelch(TurU,[],[],N,Fs,'oneside'); % turbulence spectrum
S_waveU = S_RawU-S_turU;                 % wave spectrum
[S_RawV,f] = pwelch(detrend(RawV),[],[],N,Fs,'oneside');
TurV = PhaseDec(RawV,Fs);
[S_turV,f] = pwelch(TurV,[],[],N,Fs,'oneside');
S_waveV = S_RawV-S_turV;
S_wave = S_waveU+S_waveV;
S_wave(find(S_wave <= 0)) = 10^(-16);
iw = find(f>=lf & f<=hf);
Uw = sqrt(2*sum(S_wave(iw)*(f(2)-f(1))));