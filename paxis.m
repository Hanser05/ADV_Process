function [theta,d,zr]=paxis(z)
% ----------------------------------------
% principal axis 
% [theta,d,zr]=paxis(z) 
% z = [u,v] input velocity
% zr = [u',v'] rotated velocity components
% theta = rotated angle (from x-axis)
% d = [semi-major semi-minor axis]
% ----------------------------------------
[vv,d,zr] = eof(z) ;
theta = atan2(vv(2,1),vv(1,1));
% SDATA.1