function [X,S,numIters]=NNQPSMOMulti(H,G)
% l1QPSMO for multiple columns
% Yifeng Li
% Mar. 14, 2013

[m,n]=size(G);
X=zeros(m,n);
S=zeros(m,n);
numIters=zeros(n,1);
for i=1:n
   [X(:,i),S(:,i),numIters(i)]=NNQPSMO(H,G(:,i)); 
end
end
