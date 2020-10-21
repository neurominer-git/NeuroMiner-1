% =========================================================================
% FORMAT function [param, model] = nk_GetParam_GLMFIT(Y, label, ~, ~, ModelOnly)
% =========================================================================
% Train GLMFIT models and evaluate their performance using Y & label, SlackParam & KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [param, model] = nk_GetParam2_GLMFIT(Y, label, ModelOnly, ~)
global EVALFUNC MODEFL

label(label==-1) = 0; 
warning('off')
if strcmp(MODEFL,'classification')
    model.beta  = glmfit(Y, [label ones(size(label))], 'binomial', 'link', 'logit');
else
    model.beta  = glmfit(Y, label);
end    
warning('on')
if ~ModelOnly
    Z = model.beta(1) + Y * model.beta(2:end); 
    %param.dec_values = sigmoid(model.beta(1) + Y * model.beta(2:end));
    if strcmp(MODEFL,'classification')
       param.dec_values = 1 ./ (1 + exp(-Z)); % Logistic function
    else
       param.dec_values = Z;
    end
    param.target = sign(param.dec_values);
    %param.dec_values = 1 ./ (1 + exp(-Z)); % Logistic function
    param.val = feval(EVALFUNC, label, param.target);
else
    param = [];
end
