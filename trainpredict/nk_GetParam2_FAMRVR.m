% =========================================================================
% FORMAT function [param, model] = nk_GetParam_MVTRVR(Y, label, ModelOnly, KernParam)
% =========================================================================
% Train MikRVM models and evaluate their performance using Y & label, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [param, model] = nk_GetParam2_FAMRVR(Y, label, ModelOnly, KernParam)
global SVM 

param = [];

% Generate training kernel for RVR analysis
switch SVM.kernel.kernstr
    case {' -t 0','lin','linear','lin_kernel'}
        PHI = sbl_kernelFunction(Y, Y, '+lin', 1);
    otherwise
        PHI = sbl_kernelFunction(Y, Y, ['+' SVM.kernel.kernstr ], KernParam);
end

% Learn RVR model
%switch MODEFL
    
%    case 'classification'
%        tlabel = nk_MakeRVMTargetMatrix(label);
%        [model.probs, model.weights, model.alpha, model.rv, model.kernels_] = mc_rvm(PHI, tlabel, SVM.kernel.kernstr, SVM.MVTRVR.iter);
        
%    case 'regression'
        
        % Learn RVR model
        [model.used, model.alpha, model.Mu, model.invSigma, model.Omega, model.nIters] = fmrvr(PHI, label, SVM.FAMRVR.iter, SVM.FAMRVR.tolerance);
        %[model.weights, model.rv, model.alpha, model.beta] = mvrvm(PHI, label, SVM.MVTRVR.iter);
%end

% Save model settings
model.KernParam = KernParam; model.totalSV = 0;