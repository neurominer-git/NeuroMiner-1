% =========================================================================
% FORMAT param = MAERR(expected, predicted)
% =========================================================================
% Compute Mean Absolute Error of regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2017

function [param, stdparam] = MAERR(expected, predicted)
if isempty(expected), param = []; return; end
v = predicted-expected;
param = mean(abs(v));
stdparam = std(abs(v));

end