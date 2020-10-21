% =========================================================================
% FORMAT param = SENSITIVITY(expected, predicted)
% =========================================================================
% Compute classification specificity:
% Specificity = TP / (TP + FN) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2016

function param = SENSITIVITY(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TP = sum( predicted > 0 & expected > 0 );
FN = sum( predicted < 0 & expected > 0 );

param = TP / ( TP + FN );
param = param * 100;