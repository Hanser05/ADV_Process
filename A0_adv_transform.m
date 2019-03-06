%% Transform ADV data from beam coordinate to ENU coordinate
% Reunderstand the defination of pointing down, after compared with
% ADCP&AD2CP data, we do not need to change sign if the instrument is
% pointing down!! By ltf 2018
% Ligthly modified by xjl 2016
% Based on the code developed by Geroge

clear all;clc;tic;
adv_data = load('N3181212.dat');  % Local data file
u = adv_data(:,3); v = adv_data(:,4); w = adv_data(:,5);
Vbeam = [u v w];

%% Prepare the input parameters
% Transform matrix, which can be found in *.hdr file
% T varies between different ADVs
T = [2.6982   -1.3901  -1.3076;...
    -0.0505  2.3599   -2.3184;...
    0.3433   0.3547   0.3337];

% Averaged heading, pitch and roll, the unit is degree
% Remember remove the begin and end values
heading = 173.8; pitch = -5.28; roll = 0.46;
HPR = [heading pitch roll];

%% Transform
% Call the BEAM2ENU function
Venu  = BEAM2ENU(Vbeam,HPR,T,1);

% check
figure;
plot(sqrt(u.^2+v.^2),'r');hold on
plot(sqrt(Venu(:,1).^2+Venu(:,2).^2),'b');
legend('U-BEAM','U-ENU');

%% Save in .mat
adv_data(:,3:5) = Venu;
save('adv_enu.mat','adv_data');  % If the data is big, try use -v7 matfile
% save adv4841_transformed.dat -ascii Venu
toc;




 