function [Venu,Vbeam_retransform]=BEAM2ENU(Vbeam,HPR,T,statusbit0)
%
% function Venu=Beam2ENU(Vbeam,HPR)
%
% Beam2ENU transforms Nortek ADV beam data to ENU coordinates. 
% This version assumes that the heading, pitch and roll are fixed.
% heading,pitch,roll取平均值就好，在一次观测中，这三个值（即仪器的姿态）基本上不会有太大的变化（*.sen文件）
%
% INPUTS: Vbeam = the 3D velocity in beam coordinates.  
%         HPR = [heading pitch roll]
%         heading = sensor heading in degs
%         pitch   = sensor pitch in degs
%         roll    = senor roll in degs
%
% T is the transformation matrix for beam to xyz coordinates,(in .hdr file) 
% statusbit0: if instrument is pointing down, statusbit0 equals to 1
% OUTPUTs:
%        Venu: velocity in enu coordinate
%        Vbeam_retransform: to check the validity of the transformation

% Code based on a function provided by Lee Gordon, NortekUSA(Transform.m) and George
% Voulgaris, USC(BEAM2ENUv1.m)
% Slightly modified by JLX
% Check the defination of statusbit0, slightly modified by ltf

deg2rad=pi/180;
Venu=[];
HPR=HPR*deg2rad;   %Convert to rads
%
% If instrument is pointing down (bit 0 in status equal to 1)
% rows 2 and 3 must change sign
% NOTE: For the Vector the instrument is defined to be in 
%       the UP orientation when the communication cable 
%       end of the canister is on the top of the canister, ie when 
%       the probe is pointing down.
%       For the other instruments that are used vertically, the UP 
%       orientation is defined to be when the head is on the upper
%       side of the canister.
% After test, the step below should not be done.

% If pointing up (the communication cable is lower, the probe is higher)
% change the sign of 2nd and 3th beam
if statusbit0 == 0
   T(2,:) = -T(2,:);   
   T(3,:) = -T(3,:);   
end
 

%
% Put data in columns for the calculation
[nr,nc] = size(Vbeam);      % find size of Vbeam
car=1;
if(nr>nc),
  Vbeam=Vbeam';
  car=0;
  end;
%
Vxyz = T*Vbeam;  % Beam coordinates ro XYZ coordinates
% 
     hh=HPR(1)-pi/2; 
     pp=HPR(2);
     rr=HPR(3);
%
% Make heading matrix
     H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
%
% Make tilt matrix
%
     P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
           0             cos(rr)          -sin(rr);  ...
         sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
%
% Make resulting transformation matrix
%
     R = H*P; 
%     
Venu = R*Vxyz;
%
% Notes the same program can be used to invert from ENU to XYZ and then
% Beam coordinates if use the following:
%
% Vxyz = inv(R)*Venu
% Vbeam = inv(T)*Vxyz
Vbeam_retransform = inv(R*T)*Venu;
Vbeam_retransform = Vbeam_retransform';
% return array to column/row format of Vbeam
if (car==0), Venu=Venu'; end;
end

% Principles:
% RR = H*P*T
% ENU = RR*Beam           Beam → ENU
% Beam = inv(RR)*ENU      ENU → Beam
% XYZ = T*Beam            Beam → XYZ
% Beam = inv(T)*XYZ       XYZ → Beam
% XYZ = T*inv(RR)*ENU     ENU → XYZ
