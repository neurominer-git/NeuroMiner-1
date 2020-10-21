function [rs, ds] = nk_GetTestPerf_FAMRVR(X, tXtest, ~, md, Features, k)

global SVM 

% Get Training Data
tX = nk_ExtractFeatures(X, Features, [], k);

% Generate test data kernel for RVM analysis
switch SVM.kernel.kernstr
    case {'lin','linear',' -t 0', 'lin_kernel'}
        PHI = sbl_kernelFunction(tX, tXtest, '+lin', 1);
    otherwise
        PHI = sbl_kernelFunction(tXtest, tX, ['+' SVM.kernel.kernstr ], md.KernParam);
end
% Compute the inferred prediction function
rs = PHI(:,md.used)*md.Mu;

ds = rs;

end