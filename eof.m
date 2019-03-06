function [vv,d,z]=eof(uu)
%EOF  Empirical orthogonal functions.
% [v,d,z]=eof(u)
% time series components stored as column vectors in u.
% v is the matrix of eigenvectors  特征向量
% d is a diagonal matrix of eigenvalues 特征值
% z is the time series of mode amplitudes, stored as column vectors
%
% Reference: "Wallace and Dickinson, 1972, Journal of Applied
% Meteorology, V 11, N 6, p 887-892.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.0 (12/4/96) Rich Signell (rsignell@usgs.gov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ug=denan(u);  	      %get rid of rows with missing data 
[vv,d] = eig(cov(uu));    %eigenvectors and eigenvalues of covariance matrix 协方差矩阵
d = diag(d,0);        
[d,i] = sort(d);        %sorts eigenvalues into ascending order
i = i(length(d):-1:1);  %indices of sort
d = d(length(d):-1:1);  %change to decending order
vv = vv(:,i);             %rearrange eigenvectors to correspond to eigenvalues
z = uu*vv;                % calculate mode amplitudes
