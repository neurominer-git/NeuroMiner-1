clear

load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/smallRoundBlueCellTumorOFChildhood.mat');
dataStr='SRBCT';
D=Data;
clear('Data');

% normalize D
D=normc(D);
% CV
kfold=4;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,~indTest);
testSet=D(:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);

% option
option.lambda=0.2;
option.SRMethod='l1nnls';
option.kernel='rbf';
option.param=2^0;
k=10;

% learn dictionary
trainTrain=computeKernelMatrix(trainSet,trainSet,option);
 [AtA,Y,numIter,tElapsed,finalResidual]=KSRDL([],trainTrain,k,option);
 
% sparse coding
trainTest=computeKernelMatrix(trainSet,testSet,option);
testTest=computeKernelMatrix(testSet,testSet,option);
testTest=diag(testTest);
AtTest=Y'\trainTest;
[X,tElapsed]=KSRSC(AtA,AtTest,testTest,option);


