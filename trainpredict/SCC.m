% =========================================================================
% FORMAT param = SCC(expected, predicted)
% =========================================================================
% Compute Squared Correlation Coefficient, also called Coefficient of
% Determination or alternatively % Explained Variance
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 04/2012

function param = SCC(expected, predicted)
if isempty(expected), param = []; return; end
I = ~isnan(predicted);
%param = (1 - ( sum( (expected-predicted).^2 ) / sum( (expected-mean(expected)).^2 ) )) * 100;
param = CC(expected(I),predicted(I)); param = param ^ 2 *100;
if isnan(param), error('Prediction algorithm returned non-finite performance measure'); end

end