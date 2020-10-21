% =========================================================================
% FORMAT param = RMSD(expected, predicted)
% =========================================================================
% Compute Root of Mean Squared Deviation (RMSD) of regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = RMSD(expected, predicted)
if isempty(expected), param = []; return; end
total = size(expected,1);
param =  sqrt(sum((expected-predicted).^2)/total);

end