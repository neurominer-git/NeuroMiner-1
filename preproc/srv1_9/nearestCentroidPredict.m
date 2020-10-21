function [predicted,vals]=nearestCentroidPredict(model,testSet)
% predict
% Yifeng Li
% September 17, 2012

numTe=size(testSet,2);
Kte=zeros(numTe,1);
for i=1:numTe
   Kte(i)= computeKernelMatrix(testSet(:,i),testSet(:,i),model.option);
end
vals=zeros(numTe,model.numCl);
for i=1:model.numCl
    modelThis=model.eachClass{i};
    S=computeKernelMatrix(testSet,modelThis.train,model.option);
    vals(:,i)= Kte + modelThis.mutXtXmu - 2*S*modelThis.mu;
end
% softmax
[mVal,predicted]=min(vals,[],2);
predicted=predicted-1;
end