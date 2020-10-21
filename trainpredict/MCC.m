% =========================================================================
% FORMAT param = MCC(expected, predicted)
% =========================================================================
% Compute Matthew' correlation coefficient of classification:
% MCC = (TP * TN - FP * FN) / sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = MCC(expected, predicted)

ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

param = (TP * TN - FP * FN) / sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) );

if isnan(param), param = 0; end

end