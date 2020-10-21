clc
clear
load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/smallRoundBlueCellTumorOFChildhood.mat');
dataStr='SRBCT';
D=Data;
clear('Data');

% normalize D
D=normc(D);

% reset the random number generator
s = RandStream('swb2712','Seed',1);
RandStream.setDefaultStream(s);

% CV
kfold=3;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,~indTest);
testSet=D(:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);

% compute kernel
option.kernel='rbf'; % can be linear, rbf, sigmoid, polynomial
option.param=2^0; % the parameter of kernel function
H=computeKernelMatrix(trainSet,trainSet,option);
G=computeKernelMatrix(trainSet,testSet,option);
B=computeKernelMatrix(testSet,testSet,option);

% KSRSC sparse coding
optionKSRSC.lambda=0.1;
optionKSRSC.SCMethod='l1qpAS'; % can be nnqpAS, l1qpAS, nnqpIP, l1qpIP, l1qpPX, nnqpSMO, l1qpSMO
optionKSRSC.iter=200;
optionKSRSC.dis=0;
optionKSRSC.residual=1e-4;
optionKSRSC.tof=1e-4;
X=KSRSC(H,G,diag(B),optionKSRSC);
X



