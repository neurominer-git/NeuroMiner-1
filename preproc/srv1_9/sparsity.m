function s=sparsity(Y,tol)
% Calculate the sparsity of a matrix
% Y: matrix
% s: scalar, the sparsity of Y
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 26, 2011

if nargin<2
   tol=0.001; 
end
Y=abs(Y); % absolute values
maxy=max(max(Y));
s=sum(sum(Y<=tol*maxy))/(size(Y,1)*size(Y,2));
end