%% Step3: calculate wave direction
tic; clear all; clc;
data = load('advClean_despiking.mat');
u = data.new_data2(:,4);
v = data.new_data2(:,5);
T = data.new_data2(:,2);
burst =data.new_data2(:,1);
gap = find(diff(burst)~=0);
gap = [0;gap;length(burst)];
for i = 1:length(gap)-1
    index = gap(i)+1:gap(i+1);
    Um(i,1) = mean(u(index));
    Vm(i,1) = mean(v(index));
    Tm(i,1) = mean(T(index));
end
ave_v=Vm;  ave_u=Um;
dir = atan(ave_v./ave_u)*180/pi;
for i = 1:length(ave_u)
    if ave_u(i)>0&&ave_v(i)>0
        dir(i) = 90 - dir(i) ;
    else if ave_u(i)>0&&ave_v(i)<0
            dir(i) = 90 - dir(i);
        else if ave_u(i)<0&&ave_v(i)<0
                dir(i) = 270 - dir(i);
            else if ave_u(i)<0&&ave_v(i)>0
                    dir(i) = 270 - dir(i);
                end
            end
        end
    end
end
figure;
plot (Tm,dir,'.');datetick('keepticks');
aveUVW = ones (length(Um),3)-1
aveUVW(:,1) = Tm;
aveUVW(:,2) = Um;
aveUVW(:,3) = Vm;
aveUVW(:,4)= dir;
save('adv_despiking_aveUVW.mat','aveUVW');
