function [sampleSelected,classesSelected]=sampleSelNNLS(trainSet,trainClass,testSet,option)
% sample selection by NNLS
% Yifeng Li
% Nov 18, 2011

sampleSelected=[];
classesSelected=[];
ind=crossvalind('Kfold',trainClass,option.kfold);
for i=1:option.kfold
    indCur=(ind==i);
    trainSubset=trainSet(:,indCur);
    trainSubClass=trainClass(indCur);
    fprintf('running the %d-th nnls...\n',i);
    [~,~,Y]=nnlsClassifier(trainSubset,trainSubClass,testSet,testClass,option);
    fprintf('the %d-th nnls finished\n',i);
    Y=normcl1(Y);
    Y=(Y>=option.selectThreshold); % indices of selected samples of each test sample
    Y=any(Y,2); % indices of selected samples of all test samples
    if sum(Y)==0
        warning('no sample selected, please lower down the threshold.');
    end
    sampleSelected=[sampleSelected,trainSubset(:,Y)]; % include the selected samples
    classesSelected=[classesSelected;trainSubClass(Y)];
end
end