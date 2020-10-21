function A = nk_Diversity(E, L, m, n)
% function A = nk_Diversity(P, vec, m, n)
% Compute entropy-based measure of ensemble diversity

ind0 = L ~= 0;
E = E(ind0,:); L = L(ind0);

if ~exist('m','var'), m = size(E,1); end
if ~exist('n','var'), n = size(E,2); end

E = sign(E);

rL = repmat (L, m, n);  
cE = E == rL;           % correct predictions == 1 (otherwise == 0)
cjE = sum(cE,2)*100/n;  % percent correct predictions per observation
lcjE = cjE >= 25 & cjE < 75; 

% Compute diversity
A = sum(lcjE) / m;

return