function [model] = ml_unsupervised_dimRedPCA(X, options)
% ml_unsupervised_dimRedPCA(X,options)
%
% Description:
%	 - Applies PCA to find the best low-rank approximation to 
%      the input matrix X 
%
% Options:
%   - minVariance: the minimum amount of variance to be explained by the
%                   selected components (default: 1)
%   - maxComponents: the largest number of eigenvectors to be retained
%                    in the final basis, subject to minVariance constraint
%                   (default: nVars of input X)
%
% Authors:
% 	 - Mark Schmidt (2014)
if nargin < 2
    options = [];
end

[nInstances,nVars] = size(X);

[minVariance, maxComponents, verbose] = myProcessOptions(options, ...
                                                         'minVariance', 1,...
                                                         'maxComponents',nVars, ...
                                                         'verbose', 0);

% Standardize Columns
[X,mu,sigma] = standardizeCols(X);

% Form empirical covariance
C = (1/nInstances)*X'*X;

% Compute Eigenvalues/vectors
[U,S,V] = svd(C);

% Decide which eigenvectors to keep
explained = cumsum(diag(S))/sum(diag(S));
nComponents = min(min(find(explained >= minVariance)),maxComponents);
if verbose
    fprintf('Number of Components selected: %d\n',nComponents);
    fprintf('Variance explained by basis: %.2f\n',explained(nComponents));
end

% Form Basis
basis = V(:,1:nComponents);

model.mu = mu;
model.sigma = sigma;
model.nComponents = nComponents;
model.basis = basis;
model.reduceFunc = @reduce;
end

function [Xreduced] = reduce(model,X)
    X = standardizeCols(X,model.mu,model.sigma);
    Xreduced = X*model.basis; % Note: basis is orthogonal so inverse is transpose
end