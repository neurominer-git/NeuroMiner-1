% =========================================================================
% FORMAT param = BAC2(expected, predicted)
% =========================================================================
% Compute Balanced Accuracy of classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 11/2011

function param = BAC2(expected, predicted)
global VERBOSE
if isempty(expected), param = []; return; end

ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

sens = TP / ( TP + FN ); 
spec = TN / ( TN + FP );

if ~any(expected>0)
    if VERBOSE; fprintf('\n');cprintf('red','No positive labels: Specificity will be used instead of Sensitivity * Specificity'); end
    param = spec*100;
elseif ~any(expected<0)
    if VERBOSE; fprintf('\n');cprintf('red','No negative labels: Sensitivity will be used instead of Sensitivity * Specificity'); end
    param = sens*100;
else
    param = spec * sens * 100;
end

end