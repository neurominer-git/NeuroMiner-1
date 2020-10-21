% ==========================================================================
% FORMAT [param, model] = nk_GetParam_DECTRE(Y, label, SlackParam, ModelOnly, ...
%                                           Param)
% ==========================================================================
% Train Decision Tree models (and optionally evaluate their performance on
% the training data)
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2012

function [param, model] = nk_GetParam2_DECTRE(Y, label, ModelOnly, Param)
global EVALFUNC                          

param = [];
Param(2) = 1;

% model.perc = Param(1);
% model.binthr = percentile(Y(:,2),Param(1));

% opt.Verbose = 1;
% opt.ShowPlots = 0;
if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    
    model = fitctree(Y, label,'MinLeaf',Param(1), 'MaxNumSplits', Param(2), 'PredictorSelection', 'curvature');
    %model = fitctree(Y, label,'OptimizeHyperparameters','auto', 'HyperparameterOptimizationOptions', opt);
    %fprintf('\n%s',cmdstr)
    if ~ModelOnly
        [param.target] = predict_liblin(label, Y, model);
        param.dec_values = param.target; 
        param.val = feval(EVALFUNC, label, param.dec_values);
    end
end