function [X,numIters]=l1QPProximalMulti(AtA,AtB,BtB,option)
% proximal method for l1QP
% April 11, 2012

k=size(AtA,1); % number of dictionary atoms
p=size(AtB,2); % number of new signals
X=zeros(k,p);
numIters=zeros(p,1);
finalResiduals=zeros(p,1);
for i=1:p
   [X(:,i),numIters(i),finalResiduals(i)]=l1QPProximal(AtA,AtB(:,i),BtB(i),option); 
end
end