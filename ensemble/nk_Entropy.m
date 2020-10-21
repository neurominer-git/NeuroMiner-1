function H = nk_Entropy(P, vec, m, n)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% function H = nk_Entropy(P, vec, m, n)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute entropy-based measure of ensemble diversity
%
% Inputs:
% =======
% P     = Input matrix
% vec   = unique values in input matrix
% m     = # of rows in P
% n     = # of columns in P
%
% Outputs:
% ========
% H     = Entropy of P
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2012, NeuroMiner software package

P = sign(P);
if nargin < 4, n = size(P,2); end;
if nargin < 3, m = size(P,1); end;
if nargin < 2, vec = unique(P); else vec = unique(vec(vec~=0)); end
%onevec = ones(n,1);
try  
    H = mean(arrayfun(@(i) entropy(uint8(P(:,i))),1:n));
catch
    H=0;
end
% for j=1:numel(vec)
%     % Compute frequency of k_j (class) for sample x_i
%     pf = P==vec(j);
%     f = (pf*onevec)./n;
%     f = -f.*log(f);
%     f(isnan(f)) = 0;
%     %if ~f, continue, end
%     Kx = Kx + sum(f);
% end
% H = Kx/m;
