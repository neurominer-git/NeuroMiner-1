function [A,meanVec,stdVec]=normmean0std1(A,meanVec,stdVec)
% normalize each column to have mean 0 and std 1
% Usuage: 
% [A,meanVec,stdVec]=normmean0std1(A)
% [A,meanVec,stdVec]=normmean0std1(A,meanVec,stdVec)
% A, matrix of size m by n, with samples in rows, and features in columns
% meanVec: row vector of length n
% stdVec: row vector of length n
% Example:
%     [trainSet,trainSetMean,trainSetSTD]=normmean0std1(trainSet');
%     trainSet=trainSet';
%     testSet=normmean0std1(testSet',trainSetMean,trainSetSTD);
%     testSet=testSet';
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% Sep. 25, 2010

[numR,numC]=size(A);
if nargin==1
    meanVec=nanmean(A,1);
    stdVec=nanstd(A,0,1);
    stdVec(stdVec<=eps)=eps;
end
A=(A-ones(numR,1)*meanVec)./(ones(numR,1)*stdVec);
end
