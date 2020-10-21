% ==========================================================================
% FORMAT [param, model] = nk_GetParam_WBLCOX(Y, label, ModelOnly, Tau)
% ==========================================================================
% Train Cox proportional hazards model using the Willbur base hazard rate
% and optionally evaluate its performance. If ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2020

function [param, model] = nk_GetParam2_WBLCOX(Y, label, timevec, ModelOnly, Params)
global EVALFUNC SVM 

param = [];
model = wc_train(Y, timevec, label, SVM.WBLCOX);
model.interval = SVM.WBLCOX.interval;
model.cutoff = Params;

if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_WBLCOX(Y, model);
    param.dec_values = nk_CalibrateProbabilities(param.dec_values);
    param.val = feval(EVALFUNC, label, param.dec_values);
end
