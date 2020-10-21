% ==========================================================================
% FORMAT [param, model] = nk_GetParam_MEXELM(Y, label, SlackParam, ~, ...
%                                           ModelOnly)
% ==========================================================================
% Train LIBLINEAR models and evaluate their performance using Y & label, 
% SlackParam,
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2012

function [param, model] = nk_GetParam2_MEXELM(Y, label, ModelOnly, Params)
global EVALFUNC                          

param = [];
       
if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    label(label==-1) = 2;
    if Params(2)<1, error('Number of neurons has to be >1'); end
    [model.inW, model.bias, model.outW] = mexElmTrain( Y', label, Params(2), Params(1));
  
    %fprintf('\n%s',cmdstr)
    if ~ModelOnly
        scores = mexElmPredict( inW, bias, outW, Ts );
        param.dec_values = scores(1,:)'; param.target = sign(param.dec_values);
        param.val = feval(EVALFUNC, label, param.dec_values);
    end
end