% =========================================================================
% FORMAT param = FPR(expected, predicted)
% =========================================================================
% Compute false positive rate of classification:
% False Positive Rate = 1 - TN / (TN + FP);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = FPR(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );

param = 1 - TN / (TN + FP);
param = param * 100;
