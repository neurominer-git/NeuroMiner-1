% =========================================================================
% FORMAT param = SPECIFICITY(expected, predicted)
% =========================================================================
% Compute classification specificity:
% Specificity = TN / (TN + FP)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2016

function param = SPECIFICITY(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );
     
param = TN / ( TN + FP );
param = param * 100;