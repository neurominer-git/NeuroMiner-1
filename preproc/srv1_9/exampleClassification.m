clear
load('/home/yifeng/YifengLi/Research/ML/srv1_7/data/smallRoundBlueCellTumorOFChildhood.mat');
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

% KSRSC
optionksrsc.SCMethod='l1nnlsAS';
optionksrsc.lambda=0;
optionksrsc.predicter='knn';
optionksrsc.kernel='rbf';
optionksrsc.param=2^0;
[testClassPredictedKSRSC,sparse,Y,otherOutput]=KSRSCClassifier(trainSet,trainClass,testSet,optionksrsc);
performanceKSRSC=perform(testClassPredictedKSRSC,testClass,numel(unique(trainClass)));
performanceKSRSC

% MNNLS, sub-dictionary learning, supervised dictionary learning
numMetasample=7;
optionsubdic.metaSampleMethod='vsmf'; % svd or nmf or vsmf
optionsubdic.ks=numMetasample*ones(numel(unique(classes)),1); % Colon: 3, Leukemia2: 5, Adenoma: 4, SRBCT: 5
optionsubdic.alpha2=0;
optionsubdic.alpha1=0;
optionsubdic.lambda2=0;
optionsubdic.lambda1=0; % lambda for matrix factorization
optionsubdic.t1=true;
optionsubdic.t2=true;
optionsubdic.kernelizeAY=0;
optionsubdic.method='ksrsc';
optionsubdic.kernel='linear';
optionsubdic.SCMethod='nnlsAS';
optionsubdic.lambda=2^-3; % lambda for sparse coding
optionsubdic.predictor='ns';
[metaSample,metaClass,option]=computeMetaSample(trainSet,trainClass,optionsubdic);
[testClassPredictedMNNLS1,sparse,Y,otherOutput]=KSRSCClassifier(metaSample,metaClass,testSet,optionsubdic);
performanceMNNLS1=perform(testClassPredictedMNNLS1,testClass,numel(unique(trainClass)));
performanceMNNLS1

% Use classificationTrain and classificationPredict
numMetasample=7;
optionsubdic.normalization=1;
optionsubdic.normMethod='unitl2norm';
optionsubdic.metaSampleMethod='vsmf'; % svd or nmf or vsmf
optionsubdic.ks=numMetasample*ones(numel(unique(classes)),1); % Colon: 3, Leukemia2: 5, Adenoma: 4, SRBCT: 5
optionsubdic.ifModelSelection=false;
optionsubdic.alpha2=0;
optionsubdic.alpha1=0;
optionsubdic.lambda2=0;
optionsubdic.lambda1=0; % lambda for matrix factorization
optionsubdic.t1=true;
optionsubdic.t2=true;
optionsubdic.kernelizeAY=0;
optionsubdic.method='ksrsc';
optionsubdic.kernel='linear';
optionsubdic.SCMethod='l1lsAS';
optionsubdic.lambda=2^-3; % lambda for sparse coding
optionsubdic.predictor='ns';
method='subdic';
[model,OtherOutput]=classificationTrain(trainSet,trainClass,method,optionsubdic);
[testClassPredictedMNNLS2,OtherOutput]=classificationPredict(model,testSet,testClass);
performanceMNNLS2=perform(testClassPredictedMNNLS2,testClass,numel(unique(trainClass)));
performanceMNNLS2

% Use multiClassifier
optionksrsc.normalization=1;
optionksrsc.normMethod='unitl2norm';
optionksrsc.SCMethod='l1nnlsAS';
optionksrsc.lambda=0;
optionksrsc.predicter='knn';
optionksrsc.kernel='rbf';
optionksrsc.param=2^0;
optionksrsc.search=false;
optionksrsc.ifMissValueImpute=false;

numMetasample=7;
optionsubdic.normalization=1;
optionsubdic.normMethod='unitl2norm';
optionsubdic.metaSampleMethod='vsmf'; % svd or nmf or vsmf
optionsubdic.ks=numMetasample*ones(numel(unique(classes)),1); % Colon: 3, Leukemia2: 5, Adenoma: 4, SRBCT: 5
optionsubdic.ifModelSelection=false;
optionsubdic.alpha2=0;
optionsubdic.alpha1=0;
optionsubdic.lambda2=0;
optionsubdic.lambda1=0; % lambda for matrix factorization
optionsubdic.t1=true;
optionsubdic.t2=true;
optionsubdic.kernelizeAY=0;
optionsubdic.method='ksrsc';
optionsubdic.kernel='linear';
optionsubdic.SCMethod='l1lsAS';
optionsubdic.lambda=2^-3; % lambda for sparse coding
optionsubdic.predictor='ns';

methods={'ksrsc';'subdic'};
options={optionksrsc;optionsubdic};
[testClassPredicteds,classPerforms,conMats,tElapseds,OtherOutputs]=multiClassifiers(trainSet,trainClass,testSet,testClass,methods,options);
classPerforms

% Use cvExperiment
optionksrsc.normalization=1;
optionksrsc.normMethod='unitl2norm';
optionksrsc.SCMethod='l1nnlsAS';
optionksrsc.lambda=0;
optionksrsc.predicter='knn';
optionksrsc.kernel='rbf';
optionksrsc.param=2^0;
optionksrsc.search=false;
optionksrsc.ifMissValueImpute=false;

numMetasample=7;
optionsubdic.normalization=1;
optionsubdic.normMethod='unitl2norm';
optionsubdic.metaSampleMethod='vsmf'; % svd or nmf or vsmf
optionsubdic.ks=numMetasample*ones(numel(unique(classes)),1); % Colon: 3, Leukemia2: 5, Adenoma: 4, SRBCT: 5
optionsubdic.ifModelSelection=false;
optionsubdic.alpha2=0;
optionsubdic.alpha1=0;
optionsubdic.lambda2=0;
optionsubdic.lambda1=0; % lambda for matrix factorization
optionsubdic.t1=true;
optionsubdic.t2=true;
optionsubdic.kernelizeAY=0;
optionsubdic.method='ksrsc';
optionsubdic.kernel='linear';
optionsubdic.SCMethod='l1lsAS';
optionsubdic.lambda=2^-3; % lambda for sparse coding
optionsubdic.predictor='ns';

classMethods={'ksrsc';'subdic'};
options={optionksrsc;optionsubdic};
kfold=3;
rerun=3;
fsMethod='none';
fsOption=[];
feMethod='none';
feOption=[];
[meanMeanPerformsAllRun,stdSTDAllRun,meanconMatAllRun,meantElapsedsAllRun]=cvExperiment(D,classes,kfold,rerun,fsMethod,fsOption,feMethod,feOption,classMethods,options);
meanMeanPerformsAllRun

