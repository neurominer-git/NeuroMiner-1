% ==========================================================================
% FORMAT [param, model] = nk_GetParam_MEXELM(Y, label, SlackParam, ~, ...
%                                           ModelOnly)
% ==========================================================================
% Train LIBLINEAR models and evaluate their performance using Y & label, 
% SlackParam,
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2012

function [param, model] = nk_GetParam2_RNDFOR(Y, label, ModelOnly, Param)
global EVALFUNC CMDSTR

param = [];
       
if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    
    model = CMDSTR.RFtrain(Y, label, Param(1), Param(2) );

    %fprintf('\n%s',cmdstr)
    if ~ModelOnly
        [param.target] = predict_liblin(label, Y, model);
        param.dec_values = param.target; 
        param.val = feval(EVALFUNC, label, param.dec_values);
    end
end