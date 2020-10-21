function [rs, ds] = nk_GetTestPerf_MikRVM(X, tXtest, ~, md, Features, k)
global SVM

% Get Training Data
tX = nk_ExtractFeatures(X, Features, [], k);

% Generate test data kernel for RVM analysis
switch SVM.kernel.kernstr
    case {' -t 0', 'lin','linear'}
        K = SB1_KernelFunction(tX, tXtest, SVM.kernel.kernstr, md.KernParam);
    otherwise
        K = SB1_KernelFunction(tXtest, tX, SVM.kernel.kernstr, md.KernParam);
end
% Compute the inferred prediction function
y = K * md.inference;

% Convert the output according to the likelihood (i.e. apply link function)
switch md.likelihood.InUse
    case {1,3}
        rs = y;
    case 2
        rs = SB2_Sigmoid(y)>0.5;
end

ds = rs;

end