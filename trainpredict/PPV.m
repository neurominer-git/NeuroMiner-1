% =========================================================================
% FORMAT param = PPV(expected, predicted)
% =========================================================================
% Compute false posotive rate of classification:
% Positive Predictive Value = TP / (TP + FP);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = PPV(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );

param = TP / (TP + FP);
param = param * 100;
if isnan(param), param = 0; end