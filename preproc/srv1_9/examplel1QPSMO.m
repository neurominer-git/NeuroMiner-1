% example of testing l1QPSMO

clear

load('/home/yifeng/YifengLi/Research/ML/srv1_7/data/smallRoundBlueCellTumorOFChildhood.mat');
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

% trainSet=[trainSet,-trainSet];
%testSet=testSet(:,end);

% compute kernel
option.kernel='linear';
option.param=2^0;
H=computeKernelMatrix(trainSet,trainSet,option);
G=computeKernelMatrix(trainSet,testSet,option);
B=computeKernelMatrix(testSet,testSet,option);
lambda=0;
ts1=tic;
[X1, Pset] = l1NNQPActiveSet(H, -G,lambda);
t1=toc(ts1)

ts2=tic;
[X2,S2] = NNQPSMOMulti(H, lambda-G);
t2=toc(ts2)

ts3=tic;
[X3, U] = l1QPActiveSet(H, -G,lambda);
t3=toc(ts3)

ts4=tic;
X4=l1QPSMOMulti(H,-G,lambda);
t4=toc(ts4)

ts5=tic;
lambda=0.1;
option.lambda=lambda;
option.L=size(trainSet,2);
option.residual=1e-7;
option.tof=1e-7;
X5=l1QPProximalMulti(H,G,diag(B),option);
t5=toc(ts5)







