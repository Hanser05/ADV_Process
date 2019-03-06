function [Hm1,Ha1,Ht1] = hm(p);
%直接统计得到各种波高
%Hm:最大波高
%Ha:平均波高
%Ht:1/3最大波高,即有效波高
p=p(find(p));
len=length(p);
aa=find(p<0);bb=find(p>0);

bbb=diff(bb);
ind=find(bbb>2);
j=1;
while j<length(ind)+1
    p(bb(ind(j))+1:bb(ind(j)+1)-1)=-max(-p(bb(ind(j))+1:bb(ind(j)+1)-1));
    j=j+1;
end

aaa=diff(aa);
ind=find(aaa>2);
j=1;
while j<length(ind)+1
    p(aa(ind(j))+1:aa(ind(j)+1)-1)=max(p(aa(ind(j))+1:aa(ind(j)+1)-1));
    j=j+1;
end

j=1;i=1;
while j<len
       if p(j)*p(1)>0&&p(j+1)*p(1)<0;
           wave(i)=p(1)/abs(p(1))*(p(j)-p(j+1));
           j=j+2;i=i+1;
       else j=j+1;
       end
end

ch=length(wave);
Hm1=max(wave);
wave=sort(wave,'descend');
Ha1=mean(wave);
Ht1=mean(wave(1:fix(ch/3)));
    
