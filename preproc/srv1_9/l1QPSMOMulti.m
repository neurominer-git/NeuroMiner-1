function [X,numIters]=l1QPSMOMulti(H,G,lambda)
% l1QPSMO for multiple columns
% Yifeng Li
% Mar. 14, 2013

[m,n]=size(G);
X=zeros(m,n);
numIters=zeros(n,1);
for i=1:n
   [X(:,i),numIters(i)]=l1QPSMO(H,G(:,i),lambda); 
end
end
