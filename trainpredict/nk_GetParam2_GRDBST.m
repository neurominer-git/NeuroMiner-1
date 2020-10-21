% ==========================================================================
% FORMAT [param, model] = nk_GetParam_GRDBST(Y, label, ModelOnly, Params)
% ==========================================================================
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2017

function [param, model] = nk_GetParam2_GRDBST(Y, label, ModelOnly, Params)
global EVALFUNC 

param = [];
options = nk_GenMatLearnOptions(Params);
maxIters = uint32(options.maxIters);
options.maxTreeDepth = uint32(options.maxTreeDepth);
%options.loss = 'squaredloss';
model = SQBMatrixTrain( single(Y), label, maxIters, options );

if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_GRDBST([], Y, label, model) ;
    param.val = feval(EVALFUNC, label, param.target);
end