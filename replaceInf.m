function [replaced] = replaceInf(x)   %�滻�����е�Inf, ����Χ���ݵ�ƽ��ֵ
x(isinf(x)) = NaN;   %If x is Inf, x = NaN. If not, x = x;
a = ones(size(x)+2); % a��x����������
a(a==1) = NaN;  %��a�����е�Ԫ����ΪNaN
a(2:end-1,2:end-1) = x;   
[r,c] = find(isnan(a(2:end-1,2:end-1))); %�ҵ�a��Inf���±�

for i = 1:length(r)
 bak(i) = nanmean([a(r(i),c(i)) a(r(i),c(i)+1) a(r(i),c(i)+2) a(r(i)+1,c(i))...
     a(r(i)+1,c(i)+1) a(r(i)+1,c(i)+2) a(r(i)+2,c(i)) a(r(i)+2,c(i)+1) a(r(i)+2,c(i)+2) ]);
 a(r(i)+1,c(i)+1) = bak(i);
end;
replaced = a(2:end-1,2:end-1);