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

% KSRDL
optionKSRDL.lambda=0.1;
optionKSRDL.SCMethod='nnqpAS'; % can be nnqpAS, l1qpAS, nnqpIP, l1qpIP, l1qpPX, nnqpSMO, l1qpSMO
optionKSRDL.dicPrior='uniform'; % can be uniform, Gaussian
optionKSRDL.kernel='rbf';
optionKSRDL.param=2^0;
optionKSRDL.iter=100;
optionKSRDL.dis=0;
optionKSRDL.residual=1e-4;
optionKSRDL.tof=1e-4;

H=computeKernelMatrix(trainSet,trainSet,optionKSRDL);

% training: learn the feature space
k=4; % number of basis vectors
[AtA,Y,A,pinvY]=KSRDL(trainSet,H,k,optionKSRDL); % Y is the representation of the training set in the feature space, the dictionry A is only meaningful for linear kernel. Usually we use AtA.

% prediction: project the unknown samples in the feature space
G=computeKernelMatrix(trainSet,testSet,optionKSRDL);
B=computeKernelMatrix(testSet,testSet,optionKSRDL);
AtTestSet=pinvY'*G;
X=KSRSC(AtA,AtTestSet,diag(B),optionKSRDL); % X is the representation of the test set in the feature space

