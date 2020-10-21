% =========================================================================
% FORMAT param = MSE(expected, predicted)
% =========================================================================
% Compute Mean Standand Error of regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = MSE(expected, predicted)
if isempty(expected), param = []; return; end
total = size(expected,1);
sumerr = sum((predicted-expected).*(predicted-expected));
param = sumerr/total;

end