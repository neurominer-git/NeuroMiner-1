% =========================================================================
% FORMAT function [param, model] = nk_GetParam2_matLearn(Y, label, ModelOnly, Params)
% =========================================================================
% Use Marc Schmidt's toolbox to build various ML models 
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2017

function [param, model] = nk_GetParam2_matLRN(Y, label, ModelOnly, Params)
                                            
global EVALFUNC SVM                         

matLearn_cmd = sprintf('ml_%s_%s',SVM.matLRN.learner.framework, char(SVM.matLRN.algo));

options = nk_GenMatLearnOptions(Params);

model = feval(matLearn_cmd, Y, label, options);
%if isfield(model,'w')
%    model.SV = 1-label.*(Y*model.w) >= 0; model.totalSV = sum(model.SV);
%end
param = [];
if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_matLearn([], Y, label, model) ;
    param.val = feval(EVALFUNC, label, param.target);
end
