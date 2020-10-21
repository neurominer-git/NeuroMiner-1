% =========================================================================
% FORMAT param = CC(expected, predicted)
% =========================================================================
% Compute Correlation Coefficient
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 09/2011

function param = CC(expected, predicted)
if isempty(expected), param = []; return; end
param = corrcoef(expected,predicted);
param = param(2);
if isnan(param), 
    if numel(unique(predicted))==1
        param = 0; 
    else
        error('Prediction algorithm returned non-finite performance measure'); 
    end
end
end