% =========================================================================
% FORMAT param = TPR(expected, predicted)
% =========================================================================
% Compute false posotive rate of classification:
% Positive Predictive Value = TP / (TP + FP);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2012

function param = TPR(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);
P = expected > 0;
param = sum( predicted > 0 & P ) / sum(P) * 100;

