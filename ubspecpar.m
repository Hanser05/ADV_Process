function [ubr,Tbr,iter] = ubspecpar(hs,tp,h,specform)
% UBSPECPAR！Caclulate ubr and Tbr from hs and tp using parametric spectrum
% [ubr,Tbr,iter] = ubspecform(hs,tp,h,specform)
%
% Input:
% hs！Significant wave height (m)
% tp！Peak period (s)
% h！Water depth (m)
% specform！spectral formulation to use
% specform = 'D' for Donelan spectrum (default)
% specform = 'J' for JONSWAP spectrum
% Returns:
% ubr = representative bottom orbital velocity (m/s)
% Tbr = representative bottom wave period (s)
% iter = number of iterations if Donelan spectrum is chosen
% If iter45 the calculation did not converge
% iter = 0 for the JONSWAP spectrum
%
% Patricia Wiberg, UVa
% Last modified 9 Mar 2007
if (~exist('specform','var')),
specform = 'D';
end;
g = 9.81; dffp = 0.01;
ffp = 0.2:dffp:5; %nb = length(ffp);
nt = length(tp);
iter = zeros(nt,1);
ubr = zeros(nt,1); Tbr = ubr; Tbz = Tbr;
for i = 1:nt,
m0 = hs(i).^2./16;
fp = 1./tp(i);
f = ffp.*fp; df = dffp.*fp;
kh = qkhfs(2*pi*f,h); %(seeAppendix E)
if specform=='D',
xi = 1; tol = 1e-3; m0sfD = 0.;
while abs((m0-m0sfD)/m0 >tol),
fpbar = (m0*fp^4/(g^2*6.635e-6))^(1/.7);
gam = 1.7;
if fpbar >= 0.05,
% alpha = 0.0165.*fpbar.^0.55;
if fpbar > 0.159, gam = 6.5+2.606*log(fpbar); end;
sig = 0.08+0.0013.*(fpbar).^-3;
else,
% alpha = 0.0165.*0.05.^0.55;
sig = 0.08+0.0013.*(0.05).^-3;
end;
eterm = -((ffp-1).^2)./(2.*sig.^2);
ee = exp(eterm); t2 = gam.^ee;
t1 = -(ffp).^-4;
sffpn = (ffp.^xi./ffp.^5).*exp(t1).*t2; %dimensionless spectrum
XD = 1./sum(sffpn.*dffp); %chi
sf = m0.*XD*fp.^4.*1./(f.^4*fp).*exp(t1).*t2; %Eq. 20
m0sfD = sum(sf.*df);
iter(i) = iter(i)+1;
if iter(i)>5 
    break;
end;
end;
%m1sfD = sum(f.*sf.*df);
%m2sfD = sum(f.^2.*sf.*df);
%fmD = m1sfD./m0sfD; TmD = 1./fmD;
%fzD = sqrt(m2sfD./m0sfD); TzD = 1./fzD;
su = (2*pi.*f./sinh(kh)).^2.*sf; %Eq. 5
ubrD = sqrt(2.*sum(su*df)); %Eq.6
frD = sum(su.*f.*df)./sum(su.*df); %Eq. 9
fzD = sqrt(sum(f.^2.*su.*df)./sum(su.*df));
TbrD = 1./frD; TbzD = 1./fzD;
ubr(i) = ubrD; Tbr(i) = TbrD;
end;
if specform=='J',
gam = 3.3; xi = 0; %XJ = 3.283; %The value of gam can be changed
jn = find(ffp>1);
sig = 0.07*ones(size(ffp)); sig(jn) = 0.09;
eterm = -((ffp-1).^2)./(2.*sig.^2);
ee = exp(eterm); t2 = gam.^ee;
t1 = -1.25.*(ffp).^-4;
sffpn = (ffp.^xi./ffp.^5).*exp(t1).*t2; %dimensionless spectrum
XJ = 1./sum(sffpn.*dffp); %chi
sf = m0.*fp.^4.*XJ./f.^5.*exp(t1).*t2; %Eq. 20
%m0sfJ = sum(sf.*df);
%m1sfJ = sum(f.*sf.*df);
%m2sfJ = sum(f.^2.*sf.*df);
%fmJ = m1sfJ./m0sfJ; TmJ = 1./fmJ;
%fzJ = sqrt(m2sfJ./m0sfJ); TzJ = 1./fzJ;
su = (2*pi.*f./sinh(kh)).^2.*sf; %Eq. 5
ubrJ = sqrt(2.*sum(su*df)); %Eq. 6
frJ = sum(f.*su.*df)./sum(su.*df); %Eq. 9
fzJ = sqrt(sum(f.^2.*su.*df)./sum(su.*df));
TbrJ = 1./frJ; TbzJ = 1./fzJ;
ubr(i) = ubrJ; Tbr(i) = TbrJ;
end;
end;
