function [cv, K, P] = nk_INDEPpartition(Groups, Labels, LGOflip, OutReps)
% =========================================================================
% cv = nk_INDEPpartition(Groups, Labels, LGOflip)
% =========================================================================
% This function performs a leave-group-out/-in (LGO) split of the data in 
% order to test the generalizability of a model across different groups of 
% cases. LGOflip defines whether to perform a *leave-group-out* (=0) or
% *leave-group-in* (=1) split of the data. Leave-group-in means that a model
% will be trained just on one group of subjects and tested in another or
% multiple other groups of cases, depending on the grouping vector (Group)
% Use LGO=1 only if each of your groups is of sufficient size or define 
% one vs rest scenarios using 'Groups'.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NeuroMiner 1.0, (c) Nikolaos Koutsouleris, 08/2018

if ~exist('OutReps','var') || isempty(OutReps), OutReps = 1; end
[~, Groups_u] = nk_MakeDummyVariables(Groups,'sorted');
K = numel(Groups_u); P=1;
fprintf('\nDetected %g groups in vector', K);
trainidxs = cell(OutReps,K); testidxs = cell(OutReps,K); 
for j=1:OutReps
    for i=1:K
        if LGOflip
            I = Groups == Groups_u(i);
        else
            I = Groups ~= Groups_u(i);
        end
        trainidxs{j,i} = find(I);
        indTs = find(~I & ~isnan(Labels)); testidxs{j,i} = indTs;
    end
end
cv.TrainInd = trainidxs;
cv.TestInd = testidxs;

end