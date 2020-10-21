function [rs, ds] = nk_GetTestPerf_MVTRVR(X, tXtest, ~, md, Features, k)

global SVM MODEFL

% Get Training Data
tX = nk_ExtractFeatures(X, Features, [], k);

switch MODEFL
    case 'regression'
        % Generate test data kernel for RVM analysis
        switch SVM.kernel.kernstr
            case {'lin','linear',' -t 0', 'lin_kernel'}
                PHI = SB1_KernelFunction(tX(md.rv,:), tXtest, 'lin', 1);
            otherwise
                PHI = SB1_KernelFunction(tXtest, tX(md.rv,:), SVM.kernel.kernstr, md.KernParam);
        end
        % Compute the inferred prediction function
        rs = PHI * md.weights;
        ds = rs;
    case 'classification'
        % Generate test data kernel for RVM analysis
        switch SVM.kernel.kernstr
            case {'lin','linear',' -t 0', 'lin_kernel'}
                PHI = SB1_KernelFunction(tX, tXtest, 'lin', 1);
            otherwise
                PHI = SB1_KernelFunction(tXtest, tX, SVM.kernel.kernstr, md.KernParam);
        end
        P=2; PHI2=cell(P,1);
        for p=1:P, PHI2{p}=PHI(:,md.rv{p}); end
        ds = multinomial(PHI2,md.weights,size(tXtest,1),P); ds = ds(:,1) -0.5;
        rs=sign(ds);
end

function Y=multinomial(x,W,N,P)

Y=zeros(N,P);

sum=ones(N,1);
for p=1:P
    sum=sum+ exp(x{p}*W{p});
end

for p=1:P
    Y(:,p)=exp(x{p}*W{p})./sum;
end