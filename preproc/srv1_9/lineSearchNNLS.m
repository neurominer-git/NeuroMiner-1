function [param,maxAcc]=lineSearchNNLS(trainSet,trainClass,option)
% search single parameter sigma of rbf kernel for NNLS classifier
% trainSet, matrix, each column is a training sample
% trainClass: column vector, the class labels for the training samples
% option, struct, with fields:
% option.range, vector, the range of the exponent of sigma with base 2. The
% default is [0:8], then the range of C is 2^[0:8].
% option.kfold, scalar, the number of fold of cross-validation. The default
% is 3.
% param: scalar, the optimal gamma.
% accMax, scalar, the best accuracy.
% Yifeng Li

optionDefault.range=0:8; % exponent of 2
optionDefault.kfold=3;
optionDefault.normalization=0;
optionDefault.kernel='rbf';
optionDefault.predicter='subspace';
if nargin==2
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

% save random state
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;

numRange=numel(option.range);
accs=zeros(numRange,1);
for i=1:numRange
    defaultStream.State = savedState;
    ind=crossvalind('Kfold',trainClass,option.kfold);
    acc=0;
    for k=1:option.kfold
        indTest=(ind==k);
        trainSubSet=trainSet(:,~indTest);
        testSubSet=trainSet(:,indTest);
        trainSubClass=trainClass(~indTest);
        testSubClass=trainClass(indTest);
        
%         % normalization
%         if logical(option.normalization)
%             [trainSubSet,trainSubSetMean,trainSubSetSTD]=normmean0std1(trainSubSet');
%             trainSubSet=trainSubSet';
%             testSubSet=normmean0std1(testSubSet',trainSubSetMean,trainSubSetSTD);
%             testSubSet=testSubSet';
%         end
        
        option.param=2^(option.range(i));
        testSubClassPredicted=nnlsClassifier(trainSubSet,trainSubClass,testSubSet,testSubClass,option); % call kernel NNLS classifier
        numClass=numel(unique(trainClass));
        perf=perform(testSubClassPredicted,testSubClass,numClass);
        if numClass==2
           accInd=5;
        else
            accInd=numClass+1;
        end
        acc=acc+perf(accInd);
    end
    accs(i)=acc/option.kfold;
end
[maxAcc,maxAccInd]=max(accs);
param=2^(option.range(maxAccInd));
end
