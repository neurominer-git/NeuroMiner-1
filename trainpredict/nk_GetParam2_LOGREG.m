% =========================================================================
% FORMAT function [param, model] = nk_GetParam2_LOGREG(Y, label, ModelOnly, Param)
% =========================================================================
% Train univariate logistic regression models and evaluate their performance 
% if ModelOnly = 1, return only model
% Param has to be the scalar lambda parameter
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2015

function [param, model] = nk_GetParam2_LOGREG(Y, label, ModelOnly, Param)
global EVALFUNC

nL = numel(unique(label));
model.lambda = Param;
model.theta = LRtrain(Y, label, nL, model.lambda);

if ~ModelOnly
    [param.target, param.dec_values] = LRpredict(model.theta, Y, label);
    param.dec_values = param.dec_values(:,1);
    param.val = feval(EVALFUNC, label, param.dec_values);
else
    param = [];
end
