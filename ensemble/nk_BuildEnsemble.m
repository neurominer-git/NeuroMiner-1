function [opt_F, opt_hE, opt_E, opt_D, opt_Fcat, opt_mPred] = nk_BuildEnsemble(C, L, EnsStrat, Classes, Groups)
%
% This is the main interface function for building ensembles of predictors
%
% Inputs:
% C         = Ensemble of base learners' decisions
% L         = Labels
% F         = Feature selection mask (base learner selection mask)
% EnsStrat  = Ensemble Construction Strategy
%
%
global MULTI
if ~exist('Classes','var'), Classes=[]; end;

if MULTI.flag && MULTI.train
    funcname = ['nk_Multi' EnsStrat.OptFunc];
    [opt_hE, opt_E, opt_F, opt_Fcat, opt_D, opt_mPred] = feval(funcname, C, L, EnsStrat, Classes, Groups);
else
    funcname = ['nk_' EnsStrat.OptFunc];
    [opt_hE, opt_E, opt_F, opt_D] = feval(funcname, C, L, EnsStrat);
    opt_Fcat = []; opt_mPred = [];
end

return