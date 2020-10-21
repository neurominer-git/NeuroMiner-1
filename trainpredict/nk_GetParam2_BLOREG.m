function [param, model] = nk_GetParam2_BLOREG(Y, label, ModelOnly, Param)
global SVM VERBOSE

param = [];

if ~strcmp(SVM.kernel.kernstr,'none')
    switch SVM.kernel.kernstr
        case {'lin',' -t 0', 'linear'}
            tol = Param(1);
            KernParam = 1;
        otherwise
            tol = Param(2);
            KernParam = Param(1);
    end
    model.tol = tol;
    model.kernel.KernParam = KernParam;
    model.Y = SB1_KernelFunction(Y, Y, SVM.kernel.kernstr, model.kernel.KernParam);
else
    model.tol = Param;
    model.Y = Y;
end

%t = label;
%%[ntp, d] = size(Y);
%[t,idx] = sort(t);
%x       = x(idx,:);
%x      = [ones(ntp, 1) x];

model.alpha = blogreg(model.Y, label, model.tol);
if VERBOSE, fprintf('\tB-LOG-REG: %g features ~= 0',sum(abs(model.alpha)>0)); end

if ~ModelOnly
    param.dec_values = 1./(1 + exp(-Y*model.alpha));
    param.target = sign(param.dec_values-0.5);
end

end