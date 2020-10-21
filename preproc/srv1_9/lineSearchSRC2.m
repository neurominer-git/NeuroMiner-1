function [param,maxAcc]=lineSearchSRC2(trainSet,trainClass,option)
% search single parameter lamba that weights the l_1 norm in the optimization of SRC2
% trainSet: matrix, the training set with samples in columns and features in rows.
% trainClass: column vector of numbers or string, the class labels of the traning set.
% option, struct, with fields:
% option.range, vector, the range of the exponent of lambda with base 2. The
% default is [-4:0], then the range of C is 2^[-4:0].
% option.kfold, scalar, the number of fold of cross-validation. The default
% is 3.
% option.predicter, string, the rule to interpret the sparse code, the
% default is 'subspace'.
% param: scalar, the optimal gamma.
% accMax, scalar, the best accuracy.
% Yifeng Li

optionDefault.range=-4:0; % exponent of 2
optionDefault.kfold=3;
%optionDefault.normalization=1;
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
        
        option.lambda=2^(option.range(i));
        testSubClassPredicted = SRC2(trainSubSet, trainSubClass, testSubSet, option);
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
