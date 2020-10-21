function [rs, ds] = nk_GetTestPerf_BLOREG(X, tXtest, ~, md, Features, k)
global SVM

% Get Training Data
tX = nk_ExtractFeatures(X, Features, [], k);

% Get Training Data
if isfield(md,'kernel')
    switch SVM.kernel.kernstr
        case {' -t 0', 'lin','linear'}
            K = SB1_KernelFunction(tX, tXtest, SVM.kernel.kernstr, md.kernel.KernParam);
        otherwise
            K = SB1_KernelFunction(tXtest, tX, SVM.kernel.kernstr, md.kernel.KernParam);
    end
else
    K = tXtest;
end

ds = 1./(1 + exp(-K * md.alpha));
rs = sign(ds-0.5);


end