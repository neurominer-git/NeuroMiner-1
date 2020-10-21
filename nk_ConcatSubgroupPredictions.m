function M = nk_ConcatSubgroupPredictions

% This script assumes that you have run each subgroup model on the training & CV data (the
% subgroup) and the independent test data (the rest of the population
% stored in NM.OOCV{1}). The analysis should be in NM.analysis{1}
% The following statements should concatenate the predictions across NM
% files:

% Get the NM file paths
NMfiles = spm_select(Inf,'mat','Select NM structures with subgroup analyses');
n = size(NMfiles,1);
P = cell(1,n);

% Loop through NM files and concatenate CV and OOCV cases and predictions
for i = 1:n
    load(NMfiles)
    P(i).cases = [NM.cases; NM.OOCV{1}.cases];
    P(i).predictions = [NM.analysis{1}.GDdims{1}.Regr.mean_predictions; NM.analysis{1}.OOCV{1}.RegrResults.MeanCV2PredictedValues];
end

% We will sort all subgroup predictions according to the first NM files
M = P(1).predictions;

for i=2:n
   
    % Now sort the data and grow prediction matrix by adding subgroup
    % predictions columns to the matrix
    Mi = nk_MatchID(P(i).cases,P(i).predictions,P(1).cases);
    M = [M Mi];
    
end
