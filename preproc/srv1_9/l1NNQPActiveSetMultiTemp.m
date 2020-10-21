function [X, numIters] = l1NNQPActiveSetMultiTemp(AtA, AtB,lambda)


k=size(AtA,1); % number of dictionary atoms
p=size(AtB,2); % number of new signals
X=zeros(k,p);
numIters=zeros(p,1);

for i=1:p
   [X(:,i),~,numIters(i)]=l1NNQPActiveSet(AtA,AtB(:,i),lambda); 
end
end