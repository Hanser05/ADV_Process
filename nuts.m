clear all;
clc;
DATA1 = load('CT01_02_despikingF_major_minor.mat');
DATA2 = load('CT01_01_despiking_major_minor.mat');
DATA3 = load('CT01_01_qualified_despiking.mat');
DATA4 = load('CT01_02F_qualified_despiking.mat');
DATA1.new_data1(:,1)=DATA1.new_data1(:,1)+564;
DATA4.new_data2(:,1)=DATA4.new_data2(:,1)+564;
new_data1= [DATA2.new_data1;DATA1.new_data1];
new_data2= [DATA3.new_data2;DATA4.new_data2];

save('CT01_despiking_major_minor.mat','new_data1');
save('CT01_qualified_despiking.mat','new_data2');