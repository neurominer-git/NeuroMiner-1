% =========================================================================
% FORMAT param = NETBEN(expected, predicted, pt)
% =========================================================================
% Compute net benefit based on probability threshold 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2019

function param = NETBEN(expected, predicted, pt)

if isempty(expected), param = NaN; return; end
ind0 = expected ~=0 & ~isnan(expected) & ~isnan(predicted) & predicted~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );

param = ( TP / numel(expected) ) - ( pt / (1-pt) ) * ( FP / numel(expected) ) ;