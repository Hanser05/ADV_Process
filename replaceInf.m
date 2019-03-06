function [replaced] = replaceInf(x)   %替换矩阵中的Inf, 用周围数据的平均值
x(isinf(x)) = NaN;   %If x is Inf, x = NaN. If not, x = x;
a = ones(size(x)+2); % a比x多两行两列
a(a==1) = NaN;  %将a中所有的元素设为NaN
a(2:end-1,2:end-1) = x;   
[r,c] = find(isnan(a(2:end-1,2:end-1))); %找到a中Inf的下标

for i = 1:length(r)
 bak(i) = nanmean([a(r(i),c(i)) a(r(i),c(i)+1) a(r(i),c(i)+2) a(r(i)+1,c(i))...
     a(r(i)+1,c(i)+1) a(r(i)+1,c(i)+2) a(r(i)+2,c(i)) a(r(i)+2,c(i)+1) a(r(i)+2,c(i)+2) ]);
 a(r(i)+1,c(i)+1) = bak(i);
end;
replaced = a(2:end-1,2:end-1);