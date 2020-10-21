% This is an example of how to use the NNLS, BNNLS, SRC1, SRC2, MSRC, and LRC classifiers.
clear
%% classify undercomplete data
% load data
% suppose the current folder is the one containing the SR Toolbox
load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/smallRoundBlueCellTumorOFChildhood.mat','Data','classes');
D=Data;
clear('Data');
kfold=3;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,~indTest);
testSet=D(:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);
% normalization, note: normalization is not necessary for nnls classifer
% with linear kernel
trainSet=normc(trainSet);
testSet=normc(testSet);

% classification
disp('NNLS classifier...');
tic;
[testClassPredicted,sparsity,~,otherOutput]=nnlsClassifier(trainSet,trainClass,testSet,testClass);
tElapsedNNLS=toc;
residualNNLS=otherOutput; % regression residuals
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformNNLS=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of NNLS classifier is: \n\r');
classPerformNNLS

disp('SRC2 classifier...');
tic;
% model selection
% [param,maxAcc]=lineSearchSRC2(trainSet,trainClass);
% option.lambda=param;
optionSRC2.lambda=0.001;
[testClassPredicted,sparsity,otherOutput]=SRC2(trainSet, trainClass, testSet, optionSRC2);
tElapsedSRC2=toc;
residualSRC2=otherOutput; % regression residuals
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformSRC2=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of SRC2 classifier is: \n\r');
classPerformSRC2

% classification
disp('MSRC classifier using SVD to compute meta-samples...');
tic;
unikClass=unique(trainClass);
numUnikClass=numel(unikClass);
optionMSRC.ks=5*ones(numUnikClass,1); % you can change the number of meta-sample by yourself.
optionMSRC.metaSampleMethod='svd';
[trainSetMeta,trainClassMeta]=computeMetaSample(trainSet,trainClass,optionMSRC);
size(trainSetMeta)
% model selection
% [param,maxAcc]=lineSearchSRC2(trainSet,trainClass);
% option.lambda=param;
optionMSRC.lambda=0.1;
[testClassPredicted,sparsity,otherOutput]=SRC2(trainSetMeta,trainClassMeta,testSet,optionMSRC);
tElapsedMSRCSVD=toc;
residualMSRCSVD=otherOutput; % regression residuals
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformMSRCSVD=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of MSRC classifier is: \n\r');
classPerformMSRCSVD

% classification, you may need to install my NMF toolbox (https://sites.google.com/site/nmftool/) for this
% functionality. 
disp('MNNLS classifier using NMF to compute meta-samples...');
tic;
unikClass=unique(trainClass);
numUnikClass=numel(unikClass);
optionMNNLS.ks=5*ones(numUnikClass,1); % you can change the number of meta-sample by yourself.
optionMNNLS.metaSampleMethod='vsmf'; % change 'nmf' to 'vsmf' if you do not want my NMF toolbox
[trainSetMeta,trainClassMeta]=computeMetaSample(trainSet,trainClass,optionMNNLS);
size(trainSetMeta)
[testClassPredicted,sparsity,~,otherOutput]=nnlsClassifier(trainSetMeta,trainClassMeta,testSet);
tElapsedMNNLS=toc;
residualMNNLS=otherOutput; % regression residuals
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformMNNLS=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of MSRC classifier is: \n\r');
classPerformMNNLS


disp('LRC classifier...');
tic;
[testClassPredicted,otherOutput]=lrc(trainSet,trainClass,testSet);
tElapsedLRC=toc;
residualLRC=otherOutput; % regression residuals
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformLRC=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of LRC classifier is: \n\r');
classPerformLRC

%% classify facial data
dirData='/home/yifeng/YifengLi/Research/ML/srv1_9/data/face95.mat';
load(dirData,'tensor','classes');
dataStr='face95';
D=tensor;
clear('tensor');
kfold=3;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,:,~indTest);
testSet=D(:,:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);
  
% classification
disp('NNLS classifier for facial data...');
tic;
% vectorize
numR=size(trainSet,1);
numC=size(trainSet,2);
trainSet1=matrizicing(trainSet,3);
testSet1=matrizicing(testSet,3);
trainSet1=trainSet1';
testSet1=testSet1';
% normalization, note: normalization is not necessary for nnls classifer
% with linear kernel
trainSet1=normc(trainSet1);
testSet1=normc(testSet1);
[testClassPredicted,sparsity]=nnlsClassifier(trainSet1,trainClass,testSet1,testClass);
tElapsedNNLSFace=toc;
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformNNLSFace=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of NNLS classifier is: \n\r');
classPerformNNLSFace

disp('SRC1 classifier for facial data...'); % for SRC1, you have to install L1_magic which is available at http://users.ece.gatech.edu/~justin/l1magic
% downsample
%trainSet=downsample(trainSet,16);
%testSet=downsample(testSet,16);
% vectorize
%numR=size(trainSet,1);
%numC=size(trainSet,2);
%trainSet1=matrizicing(trainSet,3);
%testSet1=matrizicing(testSet,3);
%trainSet1=trainSet1';
%testSet1=testSet1';
% normalization, note: normalization is not necessary for nnls classifer
% with linear kernel
%trainSet1=normc(trainSet1);
%testSet1=normc(testSet1);
%optionSRC1Face.epsilon=1e-3;
%optionSRC1Face.rubost=false;
%optionSRC1Face.randomCorrupt=true;
%tic;
%[testClassPredicted,sparsity,otherOutput]=src(trainSet1,trainClass,testSet1,optionSRC1Face);
%tElapsedSRC1Face=toc;
%residualsSRC=otherOutput; % regression residuals
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
%classPerformSRC1Face=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
%fprintf('The classification performance of SRC1 classifier is: \n\r');
%classPerformSRC1Face

%% classify data with missing values
% load data
% suppose the current folder is the one containing the SR Toolbox
load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/smallRoundBlueCellTumorOFChildhood.mat','Data','classes');
D=Data;
clear('Data');
kfold=3;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,~indTest);
testSet=D(:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);

% make missing values
missRate=0.7;
randNum=rand(size(trainSet));
trainSet(randNum<=missRate)=NaN;
randNum=rand(size(testSet));
testSet(randNum<=missRate)=NaN;

% classification
disp('NNLS classifier for data with missing values...');
tic;
[testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass);
tElapsedNNLSMiss=toc;
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformNNLSMiss=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of NNLS classifier is: \n\r');
classPerformNNLSMiss
%% classify overcomplete data
% suppose the current folder is the one containing the SR Toolbox
load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/breastTissue.mat','D','classes');
D=D';
classes=changeClassLabels01(classes);
kfold=4;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,~indTest);
testSet=D(:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);
% normalization, note: normalization is not necessary for nnls classifer
% with linear kernel
trainSet=normc(trainSet);
testSet=normc(testSet);

% classification
disp('NNLS classifier for overcomplete data...');
tic;
option.kernel='rbf';
% model selection
%[param,maxAcc]=lineSearchNNLS(trainSet,trainClass);
%option.param=param;
option.param=2^8;
[testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,option);
tElapsedNNLSOver=toc;
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformNNLSOver=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of NNLS classifier is: \n\r');
classPerformNNLSOver

disp('BNNLS classifier for overcomplete data...');
tic;
optionBNNLSOver.kernel='rbf';
optionBNNLSOver.kernelParamRandomAssign=true;
optionBNNLSOver.numRandom=20;
[testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,optionBNNLSOver);
tElapsedBNNLSOver=toc;
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformBNNLSOver=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of BNNLS classifier is: \n\r');
classPerformBNNLSOver

disp('SRC1 classifier...'); % for SRC1, you have to install L1_magic which is available at http://users.ece.gatech.edu/~justin/l1magic
%option.epsilon=0.05;
%option.randomCorrupt=true;
%tic;
%[testClassPredicted,sparsity]=src(trainSet,trainClass,testSet);
%tElapsedSRC1=toc;
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
%classPerformSRC1=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
%fprintf('The classification performance of SRC1 classifier is: \n\r');
%classPerformSRC1

%% classify tensor data
load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/InfBetaTensorComplete.mat','D','classes');
kfold=9;
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==1);
trainSet=D(:,:,~indTest);
testSet=D(:,:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);
% matricize the data so as to pass the valid parameters to nnls classifier
numR=size(trainSet,1);
numC=size(trainSet,2);
trainSet=matrizicing(trainSet,3);
testSet=matrizicing(testSet,3);
trainSet=trainSet';
testSet=testSet';
% normalize each sample to mean 0 and std 1
[trainSet,trainSetMean,trainSetSTD]=normmean0std1(trainSet');
trainSet=trainSet';
testSet=normmean0std1(testSet',trainSetMean,trainSetSTD);
testSet=testSet';

% classification
disp('NNLS classifier for tensor data...');
option.kernel='ds'; % dymatic systems kernel
rank=1;
lambda=5;
optionNNLSTensor.param=[numR;numC;rank;lambda];
[testClassPredicted,sparsity]=nnlsClassifier(trainSet,trainClass,testSet,testClass,optionNNLSTensor);
tElapsedNNLSTensor=toc;
% classPerformNNLS includes the accuracies of class 0, class 1, class 2, and total accuracy, and balanced accuracy
classPerformNNLSTensor=perform(testClassPredicted,testClass,numel(unique([trainClass;testClass])));
fprintf('The classification performance of NNLS classifier is: \n\r');
classPerformNNLSTensor


