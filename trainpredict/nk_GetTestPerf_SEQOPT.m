% =========================================================================
% FORMAT function [rs, ds, model] = nk_GetTestPerf_SEQOPT(Y, label, model)
% =========================================================================
% test sequence predictor model in new data 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 12/2019

function [rs, ds, model] = nk_GetTestPerf_SEQOPT(Y, L, model)

global MODEFL

[ model, ds ] = nk_OptPredSeq(Y, L, [], model);

switch MODEFL
    case 'classification'
        rs = sign(ds);
    case 'regression'
        rs = ds;
end
