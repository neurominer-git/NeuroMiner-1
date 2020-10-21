% =========================================================================
% FORMAT function [param, model] = nk_GetParam_SEQOPT(Y, label, ModelOnly, 
%                                                                 ...cmdstr)
% =========================================================================
% Train sequential optimization models and evaluate their performance using 
% Y & label, ModelOnly, cmdstr
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2020

function [param, model] = nk_GetParam2_SEQOPT(Y, label, FeatGroups, ModelOnly, PQ)

global EVALFUNC MODEFL SVM

param =[];
model = nk_OptPredSeq(Y, label, FeatGroups, [], SVM.SEQOPT.C(PQ.val(1),:), PQ.val(2), [PQ.val(3) PQ.val(4)], EVALFUNC);

if ~ModelOnly
    param.dec_values = model.optD;
    switch MODEFL
        case 'classification'
            param.target = sign(model.optD);
        case 'regression'
            param.target = model.optD;
    end
    param.val = feval(EVALFUNC, label, param.target);
end
