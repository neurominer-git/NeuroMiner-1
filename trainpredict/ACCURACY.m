% =========================================================================
% FORMAT param = ACCURACY(expected, predicted)
% =========================================================================
% Compute classification accuracy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = ACCURACY(expected, predicted)
if isempty(expected), param = []; return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);

param = sum( expected == sign(predicted) ) / numel(expected) * 100;

end