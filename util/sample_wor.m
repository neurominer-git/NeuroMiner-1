function samp=sample_wor(X,N,numsamples)
% function samp=sample_wor(X,N,numsamples)
% draws a sample WITHOUT replacement from the vector X
% each column of samp is a new random sample of size N 
% if numsamples is not not given, it is assumed to be one
if(nargin<3)
    numsamples=1;
end;
if (N>length(X))
    error('samplesize is bigger than population size\n');
end;
% make into a column vector
X=reshape(X,prod(size(X)),1);
for i=1:numsamples
    indx=randperm(length(X))';
    samp(:,i)=X(indx(1:N));
end;
