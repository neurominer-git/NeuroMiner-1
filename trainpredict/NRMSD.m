% =========================================================================
% FORMAT param = NRMSD(expected, predicted)
% =========================================================================
% Compute Normalized Root of Mean Squared Deviation (NRMSD) of regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nrmsd = 100 * [ rmse(sim, obs) / ( max(obs, na.rm = TRUE) - min(obs, na.rm = TRUE) ) ] 
% (c) Nikolaos Koutsouleris, 06/2011
function param = NRMSD(expected, predicted)
if isempty(expected), param = []; return; end
rmse = sqrt(MSE(expected,predicted));
param = (rmse / (max(expected) - min(expected))) * 100;

end