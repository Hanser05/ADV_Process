tic;
clc;clear all;close all;
%
DATA = load('GT01_02_qualified.mat');
%data=DATA.data(1:170240,:);
data=DATA.data(170273:end,:);

burst = data(:,1); 
u = data(:,4); v = data(:,5); w= data(:,6);
tol=10^-6;
%%
% u(find(abs(u)>0.4),:)=0;
% v(find(abs(v)>0.3),:)=0;
% w(find(abs(w)>0.1),:)=0;
% data(:,4)=u ;data(:,5) =v ;  data(:,6)=w;
% save('CT01_02_qualifiedplus.mat','data');

%%
% new_data2 = data;
% new_data2(:,4) = u_despiking;
% new_data2(:,5) = v_despiking;
% new_data2(:,6) = w_despiking;
% 
% 
% %%
% gordan1=find(abs(u_despiking)>abs((1000*mean(u_despiking))));
% u_despiking(:,gordan1)= nan;
% u_despiking(:,gordan1)=(u_despiking(:,gordan1+10)+u_despiking(:,gordan1-10))/2;
% gordan2=find(abs(v_despiking)>abs((1000*mean(v_despiking))));
% v_despiking(:,gordan2)=(v_despiking(:,gordan2+10)+v_despiking(:,gordan2-10))/2;
% gordan3=find(abs(w_despiking)>abs((1000*mean(w_despiking))));
% w_despiking(:,gordan3)=(w_despiking(:,gordan3+1)+w_despiking(:,gordan3-1))/2;
% %%
% gordan1=find(abs(u_despiking)>abs((1000*mean(u_despiking))));
% u_despiking(:,gordan1)= nan;
% u_despiking(:,gordan1)=(u_despiking(:,gordan1+10)+u_despiking(:,gordan1-10))/2;
% gordan2=find(abs(v_despiking)>1);
% v_despiking(:,gordan2)= nan;
% v_despiking(:,gordan2)=(v_despiking(:,gordan2+20)+v_despiking(:,gordan2-20))/2;
% gordan3=find(abs(w_despiking)>abs((1000*mean(w_despiking))));
% w_despiking(:,gordan3)=(w_despiking(:,gordan3+1)+w_despiking(:,gordan3-1))/2;


%%
new_data1 = data;
new_data1(:,4) = u_major;
new_data1(:,5) = v_minor;
new_data1(:,6) = w_despiking;
save('GT01_despiking_major_minor3.mat','new_data1');
new_data2 = data;
new_data2(:,4) = u_despiking;
new_data2(:,5) = v_despiking;
new_data2(:,6) = w_despiking;
save('GT01_qualified_despiking3.mat','new_data2');

