function [gammaOptimal,COptimal,accMax]=gridSearch(trainSet,trainClass,option)
% Grid Search for the parameters of RBF-SVM (libSVM).
% trainSet: matrix, the training set with samples in columns and features in rows.
% trainClass: column vector of numbers or string, the class labels of the traning set.
% option, struct, with fields:
% option.CRange, vector, the range of the exponent of C with base 2. The
% default is [-4:6], then the range of C is 2^[-4:6].
% option.gammaRange,the range of the exponent of gamma with base 2. The
% default is [-9:2], then the range of gamma is 2^[-9:2].
% option.kfold, scalar, the number of fold of cross-validation. The default
% is 3.
% gammaOptimal: scalar, the optimal gamma.
% COptimal: scalar, the optimal C.
% accMax, scalar, the best accuracy
% Yifeng Li
% September 05, 2011

optionDefault.CRange=-4:6; % exponent of 2
optionDefault.gammaRange=-9:2; % exponent of 2
optionDefault.kfold=3;
% optionDefault.normalization=0;

if nargin==2
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

% save random state
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;
numCRange=numel(option.CRange);
numRange=numel(option.gammaRange);
AllAccs=zeros(numCRange,numRange);
for c=1:numCRange
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
            
%             % normalization
%             if logical(option.normalization)
%                 [trainSubSet,trainSubSetMean,trainSubSetSTD]=normmean0std1(trainSubSet');
%                 trainSubSet=trainSubSet';
%                 testSubSet=normmean0std1(testSubSet',trainSubSetMean,trainSubSetSTD);
%                 testSubSet=testSubSet';
%             end
            
            C=2^(option.CRange(c));
            param=2^(option.gammaRange(i));
            
        option.trainSetting=['-s 0 -t 2 -g ',num2str(param),' -c ',num2str(C),' -v ',num2str(option.kfold)];
        option=mergeOption(option,optionDefault);
        perfAcc=svmtrain(trainSubClass,trainSubSet',option.trainSetting);
            acc=acc+perfAcc;
        end
        accs(i)=acc/option.kfold;
    end
    AllAccs(c,:)=accs;
end
[AllAccsRow,indsRow]=max(AllAccs,[],1);
[accMax,indCol]=max(AllAccsRow);
indRow=indsRow(indCol);
gammaOptimal=2^(option.gammaRange(indCol));
COptimal=2^(option.CRange(indRow));
end