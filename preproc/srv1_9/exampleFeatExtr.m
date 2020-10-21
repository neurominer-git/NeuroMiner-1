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

% feature extraction
feMethod='ksrdl';
feOption.facts=4;
optionKSRDL.lambda=0.1;
optionKSRDL.SCMethod='nnqpAS'; % can be nnqpAS, l1qpAS, nnqpIP, l1qpIP, l1qpPX, nnqpSMO, l1qpSMO
optionKSRDL.dicPrior='uniform'; % can be uniform, Gaussian
optionKSRDL.kernel='rbf';
optionKSRDL.param=2^0;
optionKSRDL.iter=100;
optionKSRDL.dis=0;
optionKSRDL.residual=1e-4;
optionKSRDL.tof=1e-4;
feOption.option=optionKSRDL;

[trainExtr,outTrain]=featureExtractionTrain(trainSet,trainClass,feMethod,feOption);
[testExtr,outTest]=featureExtrationTest(testSet,outTrain);

% classification
optionKSRSC.SCMethod='nnqpAS';
optionKSRSC.lambda=0;
optionKSRSC.predicter='knn';
optionKSRSC.kernel='rbf';
optionKSRSC.param=2^0;
[testClassPredictedKSRSC,sparse,Y,otherOutput]=KSRSCClassifier(trainSet,trainClass,testSet,optionKSRSC);
performanceKSRSC=perform(testClassPredictedKSRSC,testClass,numel(unique(trainClass)));
performanceKSRSC



