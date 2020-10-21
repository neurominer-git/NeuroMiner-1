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

% NNQP using active-set algorithm
lambda=2^-4;
XNNQPAS = NNQPActiveSet(H, lambda-G);

% l1QP using active-set algorithm
lambda=2^-3;
Xl1QPAS = l1QPActiveSet(H, -G,lambda);

% NNQP using interior-point method
option.lambda=2^-3;
XNNQPIP=NNQPIPMulti(H,-lambda+G,diag(B),option);

% l1QP using interior-point method
option.lambda=2^-3;
Xl1QPIP=l1QPIPMulti(H,G,diag(B),option);

% l1QP using proximal method
option.lambda=2^-3;
Xl1QPPX=l1QPProximalMulti(H,G,diag(B),option);

% NNQP using decomposition method/SMO
lambda=2^-4;
XNNQPSMO = NNQPSMOMulti(H, lambda-G);

% l1QP using decomposition method/SMO
lambda=2^-3;
Xl1QPSMO=l1QPSMOMulti(H,-G,lambda);








