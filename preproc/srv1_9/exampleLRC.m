clc
clear

% note: data not available now in the package
load('/home/yifeng/YifengLi/Research/ML/srv1_9/data/ionosphere.mat','ionosphere','classes'); 
dataStr='ionosphere';
D=ionosphere;
classes=changeClassLabels01(classes);

D=normc(D);
numCl=numel(unique(classes));
kfold=3;
classPerformkfold=zeros(kfold,numCl+2);
classPerformSMOkfold=zeros(kfold,numCl+2);
classPerformQPkfold=zeros(kfold,numCl+2);
for i=1:kfold
ind=crossvalind('Kfold',classes,kfold);
indTest=(ind==i);
trainSet=D(:,~indTest);
testSet=D(:,indTest);
trainClass=classes(~indTest);
testClass=classes(indTest);

% classification
disp('LRC classifier...');
numFeat=size(trainSet,1);
paramDefault=log2(sqrt(numFeat))-0;
option.kernel='rbf';
option.param=2^paramDefault;
[testClassPredicted,otherOutput]=lrc(trainSet,trainClass,testSet,testClass,option);
classPerform=perform(testClassPredicted,testClass,numCl);
classPerformkfold(i,:)=classPerform;

end
mean(classPerformkfold,1)



