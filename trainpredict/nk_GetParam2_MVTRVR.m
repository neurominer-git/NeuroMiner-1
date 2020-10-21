% =========================================================================
% FORMAT function [param, model] = nk_GetParam_MVTRVR(Y, label, ModelOnly, KernParam)
% =========================================================================
% Train MikRVM models and evaluate their performance using Y & label, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [param, model] = nk_GetParam2_MVTRVR(Y, label, ModelOnly, KernParam)
global SVM EVALFUNC MODEFL

param = [];

% Generate training kernel for RVR analysis
switch SVM.kernel.kernstr
    case {' -t 0','lin','linear','lin_kernel'}
        PHI = SB1_KernelFunction(Y, Y, 'lin', 1);
    otherwise
        PHI = SB1_KernelFunction(Y, Y, SVM.kernel.kernstr , KernParam);
end

% Learn RVR model
switch MODEFL
    
    case 'classification'
        tlabel = nk_MakeRVMTargetMatrix(label);
        [model.probs, model.weights, model.alpha, model.rv, model.kernels_] = mc_rvm(PHI, tlabel, SVM.kernel.kernstr, SVM.MVTRVR.iter);
        % Save model settings
        model.KernParam = KernParam; model.totalSV = numel(unique(cell2mat(model.rv)));
    case 'regression'
        
        % Learn RVR model
        [model.weights, model.rv, model.alpha, model.beta] = mvrvm(PHI, label, SVM.MVTRVR.iter);
        % Save model settings
        model.KernParam = KernParam; model.totalSV = length(model.rv);
end

if ~ModelOnly
    
    switch SVM.kernel.kernstr
        case {'lin','linear',' -t 0', 'lin_kernel'}
            PHI = SB1_KernelFunction(Y(model.rv,:), Y, SVM.kernel.kernstr, 1);
        otherwise
            PHI = SB1_KernelFunction(Y, Y(model.rv,:), 'rbf', KernParam);
    end
    param.target	= PHI * model.weights;
    % Compute cost function 
    param.val = feval(EVALFUNC, label, param.target);
    param.dec_values = param.target;
end

