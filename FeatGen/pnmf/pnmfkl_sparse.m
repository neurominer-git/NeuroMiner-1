function W = pnmfkl_sparse(X, r, max_iter, tol)
%
% compute PNMF based on I-divergence (non-nomarlized KL-divergence)
% input:
%   X          nonnegative data input (m times n, sparse)
%   r          number of PNMF components
%   max_iter   maximum number of iterations (defaut 5000)
%   tol        convergence tolerance (default 1e-5)
% output:
%   W          the factorizing matrix (m times r)
%
% Zhirong Yang, January 4, 2010
% 
if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter = 5000;
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-5;
end
check_step = 200;

m = size(X,1);

W = rand(m,r);

Xsum = sum(X,2);
Xnz = nonzeros(X);
for iter=1:max_iter
    W_old = W;
    if mod(iter,check_step)==0
        fprintf('iter=% 5d ', iter);
    end

    Z = sp_factor_ratio(X, W, W'*X);
    W = W .* sqrt((Z*(X'*W) +X*(Z'*W)) ...
        ./ bsxfun(@plus, Xsum'*W, bsxfun(@times, Xsum, sum(W))));
    
    diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    if diffW<tol
        fprintf('converged after %d steps.\n', iter);
        break;
    end
    
    if mod(iter,check_step)==0
        fprintf('diff=%.10f, ', diffW);
        fprintf('obj=%.10f', ...
            sum(Xnz.*log(nonzeros(sp_factor_ratio(X, W, W'*X)))) - ...
            sum(Xnz) + sum(W)*W'*Xsum);
        fprintf('\n');
    end
end
