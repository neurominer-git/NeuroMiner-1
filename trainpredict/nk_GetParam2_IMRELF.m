% =========================================================================
% FORMAT function [param, model] = nk_GetParam_IMRELF(Y, label, SlackParam, KernParam, 
%                                                                  ModelOnly)
% =========================================================================
% Train IMRELF models and evaluate their performance using Y & label, SlackParam & KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [param, model] = nk_GetParam2_IMRELF(Y, label, ModelOnly, Params)
global SVM EVALFUNC 

Y = Y';
model = SVM.IMRELF;
model.targets = label; model.X = Y;
KernParam = Params(2); SlackParam = Params(1);
if ischar(KernParam), KernParam = str2double(KernParam); end;  model.sigma = KernParam; 
if ischar(SlackParam), SlackParam  = str2double(SlackParam); end; model.lambda = SlackParam;
model.Weight = IMRelief_Sigmoid_FastImple(model.X, model.targets, model.distance, model.sigma, model.lambda, model.maxiter, 0);
if ~ModelOnly
    [param.target, param.dec_values] = IMRelief_Sigmoid_FastImple_Predict2(Y, label, Y, model.Weight, SVM.imrelief.distance, model.sigma);
    param.val = feval(EVALFUNC, label, param.target);
else
    param = [];
end