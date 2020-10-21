% =========================================================================
% FORMAT param = FSCORE(expected, predicted)
% =========================================================================
% Compute F-Score of classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2012

function param = FSCORE(expected, predicted, beta)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);
if ~exist('beta','var'), beta = 1; end
beta2 = beta^2;
TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

precision = TP / (TP+FP); recall = TP / (TP+FN);
param = (1 + beta2) * precision * recall / (beta2 * (precision + recall)) * 100;
if isnan(param), param = 0; end

end