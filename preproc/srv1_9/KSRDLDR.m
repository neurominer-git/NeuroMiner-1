function [AtA,Y,numIter,tElapsed,finalResidual]=KSRDLDR(XtX,k,option)
% kernel sparse representation based on dimension reduction
% Yifeng Li
% Feb. 16, 2012

if nargin<3
   option=[]; 
end
optionDefault.SRMethod='l1ls';
option=mergeOption(option,optionDefault);

option.kernel='linear';
XtX2=computeKernelMatrix(XtX,XtX,option);
[AtA,Y,numIter,tElapsed,finalResidual]=KSRDL(XtX2,k,option);
end