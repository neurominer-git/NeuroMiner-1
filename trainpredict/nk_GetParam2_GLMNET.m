% =================================================================================
% FORMAT function [param, model] = nk_GetParam2_GLMNET(Y, label, ModelOnly, Params)
% =================================================================================
% Use GLMNET to train an elastic net model
% See http://web.stanford.edu/~hastie/glmnet_matlab/index.html
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2017

function [param, model] = nk_GetParam2_GLMNET(Y, label, ModelOnly, Params)
                                            
global EVALFUNC SVM                     

options = nk_GenMatLearnOptions(Params);
options = glmnetSet(options);
model = glmnet(Y,label, SVM.GLMNET.family, options);

param = [];
if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_GLMNET([], Y, label, model) ;
    param.val = feval(EVALFUNC, label, param.target);
end
