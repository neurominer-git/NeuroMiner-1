function model=nearestCentroidTrain(trainSet,trainClass,option)
% nearest centroid classifier
% Yifeng Li
% September 17, 2012

optionDefault.kernel='linear';
optionDefault.param=[];
option=mergeOption(option,optionDefault);

unikCl=unique(trainClass);
numCl=numel(unikCl);
model.eachClass=cell(numCl,1);
for i=1:numCl
    ind=(trainClass==unikCl(i));
    numThisCl=sum(ind);
    modelThis.train=trainSet(:,ind);
    modelThis.K=computeKernelMatrix(modelThis.train,modelThis.train,option);
    modelThis.mu=(1/numThisCl) .* ones(numThisCl,1);
    modelThis.mutXtXmu=modelThis.mu'*modelThis.K*modelThis.mu;
    modelThis.numTr=numThisCl;
    modelThis.class=unikCl(i);
    model.eachClass{i}=modelThis;
end
model.option=option;
model.numCl=numCl;
end