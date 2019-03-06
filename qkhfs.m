function kh = qkhfs(w,h)
% QKHFS—Quick iterative calculation of kh in dispersion relationship
% kh = qkhfs( w, h )
%
% Input:
% w= Angular wave frequency = 2*pi/T where T = wave period [1/s]
% h= Water depth [m]
% Returns:
% kh= wavenumber * depth [ ]
%
% Either w or h can be a vector, but not both.
% Hard-wired for MKS units.
% Orbital velocities from kh are accurate to 3e-12 !
%
% RLSoulsby (2006)‘‘Simplified calculation of wave orbital velocities’’
% HR Wallingford Report TR 155, February 2006, Eqns. 12a–14
% csherwood@usgs.gov
% Sept 10, 2006
g = 9.81;
x = w.^2*h./g;   %从这边来看，h给的应该是平均水深，不是时间序列
y = sqrt(x).*(x<1)+x.*(x>=1);
%this appalling bit of code is about 25% faster than a for loop
t = tanh(y);
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
kh = y;
return;
