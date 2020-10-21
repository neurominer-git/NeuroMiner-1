function [X,U,numIters]=l1QPMultiTemp(H,C,lambda)
% active set l1QP for multiple columns, temp file
% Yifeng Li
% May. 28, 2013

[m,n]=size(C);
X=zeros(m,n);
U=zeros(m,n);
numIters=zeros(n,1);
for i=1:n
   [X(:,i),U(:,i),numIters(i)]=l1QP(H,C(:,i),lambda); 
end
end
