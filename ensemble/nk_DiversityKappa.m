function kappa = nk_DiversityKappa(E, L, m, n)
% function A = nk_Diversity(P, vec, m, n)
% Compute entropy-based measure of ensemble diversity

if ~exist('m','var'), m = size(E,1); end
if ~exist('n','var'), n = size(E,2); end

rL = repmat ( L, 1, n);
cE = bsxfun(@eq, E , rL);
sCE = sum(cE,2);
p = mean(sCE./m);
kappa = (1/n*sum(sCE.*(n-sCE)))/(m*(n-1)*p*(1-p));


return