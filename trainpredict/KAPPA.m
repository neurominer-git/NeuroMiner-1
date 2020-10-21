function param = KAPPA(TP, TN, FP, FN, param)
% =========================================================================
% function param = KAPPA(TP, TN, FP, FN)
% =========================================================================
% Compute Cohen's Kappa based on contigency matrix 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2015

param.X = [TP FP; FN TN];
n       = sum(param.X(:)); %Sum of Matrix elements
x       = param.X./n;
r       = sum(x,2); %rows sum
s       = sum(x); %columns sum
f       = diag(ones(1,2));
Ex      = r*s;
po      = sum(sum(x.*f)); pe = sum(sum(r*s.*f)); %expected proportion for random agree
k       = (po - pe) / (1 - pe);
wbari   = r'*f; wbarj = s*f; wbar = repmat(wbari',1,2)+repmat(wbarj,2,1);
a       = Ex.*((f-wbar).^2);
var_k   = (sum(a(:))-pe^2)/(n*((1-pe)^2)); %variance
z       = k/sqrt(var_k); %normalized kappa
p       = (1-0.5*erfc(-abs(z)/realsqrt(2)))*2;
param.Kappa = k;
param.Kappa_Z = z;
param.Kappa_P = p;