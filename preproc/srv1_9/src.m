function [testClassPredicted,sparsity,otherOutput]=src(trainSet,trainClass,testSet,option)
% sparse representation classification approach (SRC1)
% trainSet, matrix, each column is a training sample
% trainClass: column vector, the class labels for the training samples
% testSet: matrix, each column is a new or testing sample
% option, struct, with fields:
% option.epsilon, scalar, the tolerance in the stable version, the default is 0.05.
% option.predicter, string, the rule to interpret the sparse code, it can be
% 'subspace' (default),'max','kvote'.
% option.randomCorrupt, logical, if use the version for corrution/noise,
% the default is false.
% option.rubost, if run the stable/robust version of SCR1, the default is
% false.
% Yifeng Li
% August 04, 2010
% note: each sample has to be normalized to unit l2 norm before runing it

% % normalization to length 1
% trainSet=normc(trainSet);
% testSet=normc(testSet);

%optionDefault.p=16;
optionDefault.epsilon=0.05;
optionDefault.predicter='subspace';
optionDefault.randomCorrupt=false;
optionDefault.rubost=false;
if nargin<4
   option=[]; 
end
option=mergeOption(option,optionDefault);
% trainSet=downsample(trainSet,option.p);
% testSet=downsample(testSet,option.p);

if option.randomCorrupt
    trainSet=[trainSet,eye(size(trainSet,1),size(trainSet,1))];
end

% training step, obtain sparse coefficients in columns of Y
Y=zeros(size(trainSet,2),size(testSet,2));
% if option.randomCorrupt
%     Y= trainSet\testSet;
%     testSet=testSet-Y(end-size(trainSet,1)+1:end,:);
%     Y=Y(1:size(trainSet,1),:);
%     trainSet=trainSet(:,1:(size(trainSet,2)-size(trainSet,1)));
% else
%     for i=1:size(testSet,2)
%         y0=pinv(trainSet)*testSet(:,i);
%         Y(:,i)= l1qc_logbarrier(y0, trainSet, [], testSet(:,i), option.epsilon);
%     end
% end

for i=1:size(testSet,2)
    if option.randomCorrupt
        y0=pinv(trainSet)*testSet(:,i);
        yi=l1eq_pd(y0, trainSet, [], testSet(:,i),option.epsilon);
        Y(:,i)=yi;
        testSet(:,i)=testSet(:,i)-yi(end-size(trainSet,1)+1:end);
    else
        y0=pinv(trainSet)*testSet(:,i);
        if option.rubost
            Y(:,i)= l1qc_logbarrier(y0, trainSet, [], testSet(:,i), option.epsilon);
        else
            yi=l1eq_pd(y0, trainSet, [], testSet(:,i),option.epsilon);
            Y(:,i)=yi;
        end
    end
end

if option.randomCorrupt
    Y=Y(1:size(trainSet,2)-size(trainSet,1),:);
    trainSet=trainSet(:,1:(size(trainSet,2)-size(trainSet,1)));
end
% calculate sparsity
sparsity=sum(sum(abs(Y)<=0.0001))/(size(Y,1)*size(Y,2));

otherOutput=[];
% predict step
switch option.predicter
    case  'max'
        [val,ind]=max(Y,[],1);
        testClassPredicted=trainClass(ind);
    case 'kvote'
        for s=1:size(Y,2)
            [sortedCoeff,ind] = getBestScores(Y(:,s),option.k);
            predicted(s,:)=trainClass(ind);
        end
        testClassPredicted=vote(predicted);
    case 'subspace'
        [testClassPredicted,residuals]=subspace(Y,testSet,trainSet,trainClass);
        otherOutput=residuals;
end
end