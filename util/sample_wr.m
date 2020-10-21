function samp=sample_wr(X,N,numsamples)
% function samp=sample_wr(X,N,numsamples)
% draws a sample with replacement from the vector X
% each column of samp is a new random sample of size N 
% if numsamples is not not given, it is assumed to be one
if(nargin<3)
    numsamples=1;
end;
samp=X(unidrnd(length(X),N,numsamples));