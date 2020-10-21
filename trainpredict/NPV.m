% =========================================================================
% FORMAT param = NPV(expected, predicted)
% =========================================================================
% Compute Negative Predictive Value of classification: TN / (TN + FN);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2016

function param = NPV(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TN = sum( predicted < 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

param  = (TN / (TN + FN)) * 100;
if isnan(param), param = 0; end