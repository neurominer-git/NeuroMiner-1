function param = PSI(expected, predicted)

if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);
TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

PPV   = (TP / (TP + FP)) * 100;
NPV   = (TN / (TN + FN)) * 100;
param = PPV + NPV - 100;
if ~isfinite(param), param = 0; end