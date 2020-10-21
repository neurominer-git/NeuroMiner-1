% =========================================================================
% FORMAT param = BAC(expected, predicted)
% =========================================================================
% Compute Balanced Accuracy of classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = BAC(expected, predicted)
global VERBOSE

if isempty(expected), param = NaN; return; end
ind0 = expected ~=0 & ~isnan(expected) & ~isnan(predicted) & predicted~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

sens = TP / ( TP + FN ); 
spec = TN / ( TN + FP );

if ~any(expected>0)
    if VERBOSE; fprintf('\n');cprintf('red','No positive labels: Specificity will be used instead of BAC'); end
    param = spec*100;
elseif ~any(expected<0)
    if VERBOSE; fprintf('\n');cprintf('red','No negative labels: Sensitivity will be used instead of BAC'); end
    param = sens*100;
else
    param = ((spec + sens) / 2)*100;
end

end