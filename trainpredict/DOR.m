% =========================================================================
% FORMAT param = DOR(expected, predicted)
% =========================================================================
% Compute Diagnostic Odds Ratio of classification: 
% sens/(1-spec)/((1-spec)/sens);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2016

function param = DOR(expected, predicted)
if isempty(expected), param = []; return; end

sens = SPECIFICITY(expected,predicted);
spec = SENSITIVITY(expected,predicted);

param = sens/(1-spec)/((1-spec)/sens);