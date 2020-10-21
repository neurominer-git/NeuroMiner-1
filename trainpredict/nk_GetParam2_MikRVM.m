% =========================================================================
% FORMAT function [param, model] = nk_GetParam_MikRVM(Y, label, ~, KernParam, ModelOnly)
% =========================================================================
% Train MikRVM models and evaluate their performance using Y & label, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [param, model] = nk_GetParam2_MikRVM(Y, label, ModelOnly, KernParam)
global SVM EVALFUNC RVM MODEFL

param = [];
model.likelihood	= SB2_Likelihoods(RVM.LikelihoodModel);
        
% Generate kernel for RVM analysis
K = SB1_KernelFunction(Y, Y, SVM.kernel.kernstr, KernParam);

% Learn sparse Baysian model
if strcmp(MODEFL,'classification'), label(label==-1)=2; end
[P, HP, D] = SparseBayes(RVM.LikelihoodModel, K, label, RVM.UserOpt, RVM.ParamSet);

% % Manipulate the returned weights for convenience later
w_infer						= zeros(size(K,1),1);
w_infer(P.Relevant)         = P.Value;

% Compute the inferred prediction function
y = K * w_infer;

% Convert the output according to the likelihood (i.e. apply link function)
switch model.likelihood.InUse
    case {1,3}
        param.target	= y;
    case 2
        param.target    = double(SB2_Sigmoid(y)>0.5);
end

if ~ModelOnly
    % Compute cost function 
    param.val = feval(EVALFUNC, label, param.target);
    param.dec_values = param.target;
end

% Save model settings
model.P = P; model.HP = HP; model.D = D ; model.inference = w_infer;
model.KernParam = KernParam; model.totalSV = length(P.Relevant);