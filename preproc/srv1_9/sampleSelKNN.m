function [sampleSelected,classesSelected]=sampleSelKNN(trainSet,trainClass,testSet,option)

sampleSelected=[];
classesSelected=[];
inds=[];
[rTe,cTe]=size(testSet);
[rTr,cTr]=size(trainSet);
distEu=euclidean(testSet,trainSet);
for i=1:cTe
  [sortedDist,index] = getBestScores(distEu(i,:),option.neighbor);
  index=index(:);
  inds=[inds;index];  
end
inds=unique(inds);
 sampleSelected=[sampleSelected;trainSet(:,inds)];
 classesSelected=[classesSelected;trainClass(inds)];
end