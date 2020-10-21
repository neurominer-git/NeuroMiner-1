function [sl, sv, pcal, pcav, sampleMean, paths] = spca(X, Gram, options)

% SPCA The SPCA algorithm for computing sparse principal components.
% [SL SV PCAL PCAV PATHS] = SPCA(X, Gram, K, LAMBDA, STOP) computes
% sparse principal components of the data in X. X is an n x p matrix
% where n is the number of observations and p is the number of
% variables. Gram = X'X is the p x p Gram matrix. Either X, Gram or
% both may be supplied. Pass an empty matrix as argument if either of X
% or Gram is missing.
%
% Returns SL, the sparse loadings (aka principal component directions),
% SV, the variances of each component, PCAL, the loadings of a regular
% PCA, PCAV, the variances of a regular PCA and PATHS, an optional
% struct containing the loading paths for each component as functions of
% iteration number.
%
% K is the desired number of sparse principal components.
%
% LAMBDA specifies the positive ridge term coefficient. If LAMBDA is set
% to infinity, soft thresholding is used to calculate the components.
% This is appropriate when p>>n and results in a significantly faster
% algorithm.
%
% STOP is the stopping criterion. If STOP is negative, its absolute
% value corresponds to the desired number of variables. If STOP is
% positive, it corresponds to an upper bound on the L1-norm of the BETA
% coefficients. STOP = 0 results in a regular PCA. Negative STOP cannot
% be used with LAMBDA set to infinity.
%
% The extra inputs MAXITER and TRACE determine the maximum number of
% iterations and control output.
%
% Note that if X is omitted, the absolute values of SV cannot be
% trusted, however, the relative values will still be correct.
%
% Author: Karl Skoglund, IMM, DTU, kas@imm.dtu.dk
% Reference: 'Sparse Principal Component Analysis' by Hui Zou, Trevor
% Hastie and Robert Tibshirani, 2004.
%% Input checking and initialization
global VERBOSE

if isempty(VERBOSE), VERBOSE = 1; end

% Set defaults
maxiter = 300;
stop = 0;
lambda = 1e-6;
savepaths = 0;

if exist('options','var') && isstruct(options)
    if isfield(options,'maxiter'), maxiter = options.maxiter; end
    if isfield(options,'stop'), stop = options.stop; end
    if isfield(options,'lambda'), lambda = options.lambda; end
    if isfield(options,'savepaths'), savepaths = options.savepaths; end
    if isfield(options,'ReducedDim'), K = options.ReducedDim; end
end
    
if isempty(X) && isempty(Gram)
    error('Must supply a data matrix or a Gram matrix or both.');
end
if lambda == inf && any(stop < 0)
    warning('Cannot use negative STOP (number of variables criterion) with LAMBDA = inf (soft thresholding).');
end


%% SPCA algorithm setup
if isempty(X)
    [V, D] = eig(Gram);
    X = V*sqrt(abs(D))*V';
end

% Mean-centering of data
sampleMean = mean(X);
X = bsxfun(@minus,X, sampleMean);

[~, D, pcal] = svd(X, 'econ');
[n, p] = size(X);
pcav = diag(D).^2/n;
K = min([K p n-1]);
if savepaths
    for k = 1:K
        paths(k).data = [];
    end
end
A = pcal(:,1:K);
B_norm = zeros(p,K);
iter = 0;
converged = 0;
%% SPCA loop
while ~converged && iter < maxiter
    iter = iter + 1;
    if VERBOSE && ~mod(iter, 10)
        fprintf('\nIteration %g, diff = %g', iter, max(abs(B_old(:) - B_norm(:))));
    end
    B_old = B_norm;
    for j = 1:K
        if length(stop) == K
            jstop = stop(j);
        else
            jstop = stop(1);
        end
        if lambda == inf
            % Soft thresholding, calculate beta directly
            if isempty(Gram)
                AXX = (A(:,j)'*X')*X;
            else
                AXX = A(:,j)'*Gram;
            end
            b = sign(AXX).*max(0, abs(AXX) - jstop);
        else
            % Find beta by elastic net regression
            b = larsen(X, X*A(:,j), lambda, jstop, 0);
        end
        B(:,j) = b(end,:)';
    end
    % Normalize coefficients
    B_norm = sum(B.^2);
    B_norm(B_norm == 0) = 1;
    B_norm = B./sqrt(ones(p,1)*B_norm);
    converged = max(abs(B_old(:) - B_norm(:))) < 1e-6;
    % Save coefficient data
    if savepaths
        for k = 1:K
            paths(k).data = [paths(k).data B(:,k)];
        end
    end
    % Update A
    if isempty(Gram)
        [U D V] = svd(X'*(X*B), 'econ');
    else
        [U , ~, V] = svd(Gram*B, 'econ');
    end
    A = U*V';
end
%% Normalization of loadings
% Normalize coefficients such that loadings has Euclidian length 1
B_norm = sum(B.^2);
B_norm(B_norm == 0) = 1;
sl = B./sqrt(ones(p,1)*B_norm);
%% Order modes such that maximal total explained variance is achieved
ss = X*sl; % sparse scores
sv = zeros(K, 1); % adjusted variances
O = 1:K; % ordering
for k = 1:K
    ss_var = sum(ss.^2)/n; % variances of scores
    [sv(k), max_col] = max(ss_var);
    s = ss(:,max_col); % column to factor out
    s_norm = s'*s;
    if s_norm > eps,
        O(O == max_col) = O(k);
        O(k) = max_col;
        ss(:,O) = ss(:,O) - s*s'*ss(:,O)/s_norm; % factor out chosen column
    end
end
sl = sl(:,O); % change order of loadings
%% Print information

if VERBOSE
    if p < 20
        fprintf('\n\n --- Sparse loadings ---\n\n');
        disp(sl)
    end
    fprintf('\n --- Adjusted variances, Variance of regular PCA ---\n\n');
    disp([sv/sum(pcav) pcav(1:K)/sum(pcav)])
    fprintf('Total: %3.2f%% %3.2f%%\n', 100*sum(sv/sum(pcav)), 100*sum(pcav(1:K)/sum(pcav)));
    fprintf('\nNumber of nonzero loadings:\n');
    fprintf('%g',sum(abs(sl) > 0));
end