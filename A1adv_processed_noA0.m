%% Remove bad data of ADV
% Note, if your raw ADV data is already in enu coordinate(no using A0), use
% this script. If not, try A0_adv_transform first.
% Creater: CHEN Jingdong
% Lightly modified at 2019/03/06 by Hanser05

clear all;close all;clc
tic
hp = 0.246+0.9; % hp height of the pressure sensor above the bottom, unit, m  
dt = 1/16;    % dt is the sample interval in s
second = datenum([0 0 0 0 0 1]); % second: the value of one second
FILE = load('N3181212.dat');   % raw ADV data
sample = FILE(:,2); burst = FILE(:,1); p = FILE(:,15);  
u_tmp=FILE(:,3);v_tmp=FILE(:,4);w_tmp=FILE(:,5);
SNR1=FILE(:,9);SNR2=FILE(:,10);SNR3=FILE(:,11);
cor1=FILE(:,12);cor2=FILE(:,13);cor3=FILE(:,14);

%% processing time
Time = load('N3181212.vhd');   % vhd file   
sec = Time(:,6);min = Time(:,5);hor = Time(:,4);day = Time(:,2);mon = Time(:,1);yer = Time(:,3);
time = [yer mon day hor min sec]; time = datenum(time);                    
ttall = datestr(time,'mmm.dd,yyyy HH:MM:SS'); ttall = cellstr(ttall);
clear Time sec min hor day mon yer 
%% 
% build structure
% in this code, index ii was used for sample number
% use Cor and SNR to define good data 
gap = find(diff(burst)~=0);
gap = [0;gap;length(sample)];
sample_num(1:length(gap)-1,1)=0;
for ii = 1:length(gap)-1
    index = gap(ii)+1:gap(ii+1);
    sample_num(ii,1) = max(index)-min(index)+1;
    DATA(ii) = struct('time',sample(index)*dt*second+time(ii),'burst',burst(index),'pressure',p(index),...
    'u',u_tmp(index),'v',v_tmp(index),'w',w_tmp(index),...
    'cor1',cor1(index),'cor2',cor2(index),'cor3',cor3(index),...
    'SNR1',SNR1(index),'SNR2',SNR2(index),'SNR3',SNR3(index));
    a = DATA(ii).cor1>=70&DATA(ii).cor1<101; b = DATA(ii).cor2>=70&DATA(ii).cor2<101;
    c = DATA(ii).cor3>=70&DATA(ii).cor3<101;
    d = DATA(ii).SNR1>20; e = DATA(ii).SNR2>20; f = DATA(ii).SNR3>20; x = DATA(ii).pressure > 0;
    g = a&b&c&x; gg=d&e&f; ggg = g&gg; LEN_gooddata = length(find(ggg));
    if LEN_gooddata/sample_num(ii,1)>=0.9
        HEAD(ii) = struct('burstnumber',burst(index(1)),'sample_number',sample_num(ii),...
        'time',time(ii),'timestr',ttall(ii),'quality','good','tidenumber',nan,...
        'goodsamples_ratio',LEN_gooddata/sample_num(ii,1));
    else
        HEAD(ii) = struct('burstnumber',burst(index(1)),'sample_number',sample_num(ii),...
        'time',time(ii),'timestr',ttall(ii),'quality','bad','tidenumber',nan,...
        'goodsamples_ratio',LEN_gooddata/sample_num(ii,1));
    end
    %HEAD has five parameters: burstnumber, time, quality, tidenumber, and goodsamples_ratio
end
%%
clear g gg ggg a b c d e f cor1 cor2 cor3 u_tmp v_tmp w_tmp SNR1 SNR2 SNR3 second ii gap LEN_gooddata
%% check data
LEN_burst = length(HEAD);
for i = 1:LEN_burst
    Ratio(i) = HEAD(i).goodsamples_ratio;
    Pressure(i) = mean(DATA(i).pressure);
end
figure;
subplot(2,1,1);
plot(1:LEN_burst,Ratio,'o');
ylabel('Ratio of good samples');
subplot(2,1,2);
plot(1:LEN_burst,Pressure+hp,'.');

ylabel('pressure recorded (dbar)');
%% --------------------------------------------------------------------------
% find the good burst number for next step
clear index
k=1;
for i = 1:LEN_burst
    if strcmp(HEAD(i).quality,'good')
       index(k) = i;  % the number of good burst will be stored in the index
       k=k+1;
    end
end
clear k
GOODburstnumber = length(index);  
%% If ADV is on tidal flat, relatively decrease the quality control standard
% NOTE: this method need be changed if the instrument wasn't exposed to the air between two tidal cycles.
k=1; % k is used to calculate the number of tidal cyclings.
j=1;
for i = 1:GOODburstnumber-1
     if (index(i+1)-index(i))>1&&(index(i+1)-index(i))<20
         disp('there is an abnormal discontinuity point during a tidal cycles')
         disp('at burst,plz check it later!')
         disp(HEAD(index(i)).burstnumber)
         if strcmp(HEAD(index(i)).quality,'good')
             HEAD(index(i)).tidenumber = k;
             tide(j) = k;
             j = j+1;
         end
         for m = 1:index(i+1)-index(i)-1
             if strcmp(HEAD(index(i)+m).quality,'bad')&&HEAD(index(i)+m).goodsamples_ratio>0.6
                 HEAD(index(i)+m).quality='good' ;
                 HEAD(index(i)+m).tidenumber=k;
                 tide(j)=k;
                 j=j+1;
             end
        end
        continue
     end
    HEAD(index(i)).tidenumber = k;
    tide(j) = k;
    j=j+1;
    % tide: tide number for varible 'index',also used for test
    if (index(i+1)-index(i))>40      % dry bursts number lager than 40
       k = k+1;                           
    end       
end
HEAD(index(GOODburstnumber)).tidenumber = k;
tide(j) = k;
%% After above process, good burst number may increse, rewrite index
clear index
k=1;
for ii=1:LEN_burst
    if strcmp(HEAD(ii).quality,'good')
        index(k)=ii;        % the number of good burst will be stored in the index
        k=k+1;
    end
end
clear k
GOODburstnumber=length(index); 
%% Atmosphere pressure correction
j = 1;
for k = 1:max(tide)
    aa = find(tide==k);
    ii = index(aa);
    if ii(1)==1&&ii(length(aa))~= LEN_burst
        pbefore(k) = mean(DATA(ii(length(aa))+j).pressure);% calculate air pressure before the tide 
        pafter(k) = mean(DATA(ii(length(aa))+j).pressure);% calculate air pressure after the tide
        elseif ii(1)~=1&&ii(length(aa))==LEN_burst
            pbefore(k) = mean(DATA(ii(1)-j).pressure);
            pafter(k) = mean(DATA(ii(1)-j).pressure);
        elseif ii(1)~=1&&ii(length(aa))~=LEN_burst
            pbefore(k) = mean(DATA(ii(1)-j).pressure);
            pafter(k)=mean(DATA(ii(length(aa))+j).pressure);
    end
   % calculate p0
   p0(k)=0.5*(pbefore(k) + pafter(k));
end

for i = 1:GOODburstnumber
    xx = 1:sample_num(index(i));
    xx = xx(:);
    pp=polyfit(xx,DATA(index(i)).pressure,1);  
    trend=polyval(pp,xx);
   DATA(index(i)).wavepressure = DATA(index(i)).pressure - trend;
     k = HEAD(index(i)).tidenumber;      
    HEAD(index(i)).airpressure = p0(k); % acturally, this is not air pressue in this site
    HEAD(index(i)).waterdepth = mean(DATA(index(i)).pressure)+hp...
        -HEAD(index(i)).airpressure;
      DATA(index(i)).waterdepth = DATA(index(i)).pressure + hp -HEAD(index(i)).airpressure;
    % wave pressure data is stored in DATA variable
    % water depth is assigned to be the one at the MIDDLE of the
    % stored in HEAD variable
    % They only exist in "good" burst   
end
%% save data
for i=1:length(index)
    data_temp(i)=DATA(index(i));
    head_temp(i)=HEAD(index(i));
end
%save('rd02_adv_20140921-0925_data.mat','-struct','data');
%save('rd02_adv_20140921-0925_head.mat','-struct','head');
% as structure data is not an array, may not be easily plotted or saved, so
% change it into array
%%
sum=0;

for i=1:length(index)
    data(sum+1:sum+sample_num(index(i)),1)=data_temp(i).burst;
    data(sum+1:sum+sample_num(index(i)),2)=data_temp(i).time;
    data(sum+1:sum+sample_num(index(i)),3)=data_temp(i).pressure;  %这里应该存校对后的水深吧，不然后面又要校对s
    data(sum+1:sum+sample_num(index(i)),4)=data_temp(i).u;
    data(sum+1:sum+sample_num(index(i)),5)=data_temp(i).v;
    data(sum+1:sum+sample_num(index(i)),6)=data_temp(i).w;
    data(sum+1:sum+sample_num(index(i)),7)=data_temp(i).SNR1;
    data(sum+1:sum+sample_num(index(i)),8)=data_temp(i).SNR2;
    data(sum+1:sum+sample_num(index(i)),9)=data_temp(i).SNR3;
    data(sum+1:sum+sample_num(index(i)),10)=data_temp(i).wavepressure;
    data(sum+1:sum+sample_num(index(i)),11)=data_temp(i).waterdepth;
    sum=sum+sample_num(index(i));
    tttt(i,1) = mean(data_temp(i).time);
    check_p(i,1) = mean(data_temp(i).pressure);
end
for i=1:length(index)   
     head(i,1)=head_temp(i).time;           
     head(i,2)=head_temp(i).goodsamples_ratio;    %the good sample ratio of this burst
     head(i,3)=head_temp(i).waterdepth;
     head(i,4)=head_temp(i).tidenumber;
end
 wavepressure = rmfield(data_temp, { 'u', 'v', 'w', 'cor1', 'cor2', 'cor3', 'SNR1', 'SNR2', 'SNR3'});

%%

save('advClean.mat','data');
save('advClean_waterDepth.mat','head');
save('advClean_wavepressure.mat','wavepressure');

toc