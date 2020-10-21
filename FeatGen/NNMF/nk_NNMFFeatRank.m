function [trainExtr, testExtr, outTrain, outTest] = nk_NNMFFeatRank(Y, label, Ynew)

feMethod='nmf';
optionFE.facts=15;
Y=Y';

[trainExtr,outTrain]=featureExtractionTrain(Y,[],label,feMethod,optionFE);

if exist('Ynew','var') && ~isempty(Ynew) && nargout==4
    Ynew = Ynew';
    [testExtr,outTest]=featureExtrationTest(Y,Ynew,outTrain);
else
    testExtr=[]; outTest=[];
end