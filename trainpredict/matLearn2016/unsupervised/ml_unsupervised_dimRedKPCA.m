function [model] = ml_unsupervised_dimRedKPCA(X, options)
% ml_unsupervised_dimRedKPCA(X,options)
%
% Description:
%	 - Uses kernel PCA to find the best low-rank approximation to 
%      the input matrix X using a user-specified kernel
%
% Options:
%   - minVariance: the minimum amount of variance to be explained by the
%                   selected components (default: 1)
%   - maxComponents: the largest number of eigenvectors to be retained
%                    in the final basis, subject to minVariance constraint
%                   (default: nVars of input X)
%   - kernelFunc: choice of kernel representation of basis (default:
%                 linear)
%   - kernelArgs: cell array or struct containing arguments for the kernel
%
% Authors:
% 	 - Mark Schmidt (2014)

if nargin < 2
    options = [];
end

[nInstances,nVars] = size(X);

[minVariance,maxComponents,kernelFunc,kernelArgs] =...
    myProcessOptions(options,'minVariance',1,'maxComponents',nVars,...
                     'kernelFunc',@ml_kernel_gram,'kernelArgs',{});

% Standardize Columns
[X,mu,sigma] = standardizeCols(X);

% Form Gram Matrix
K = kernelFunc(X,X,kernelArgs);

% Subtract off mean in feature space
N = (1/nInstances)*ones(nInstances);
K = K - N*K - K*N + N*K*N;

% Compute Eigenvalues/vectors
[U,S,V] = svd(K);

% Decide which eigenvectors to keep
S = S/nInstances;
explained = cumsum(diag(S))/sum(diag(S));
nComponents = min(min(find(explained >= minVariance)),maxComponents);
fprintf('Number of Components selected: %d\n',nComponents);
fprintf('Variance explained by basis: %.2f\n',explained(nComponents));

% Form Basis
basis = V(:,1:nComponents);

model.mu = mu;
model.sigma = sigma;
model.nComponents = nComponents;
model.basis = basis;
model.Xtrain = X;
model.nTrain = nInstances;
model.reduceFunc = @reduce;
model.kernelFunc = kernelFunc;
model.kernelArgs = kernelArgs;
end

function [Xreduced] = reduce(model,X)
    X = standardizeCols(X,model.mu,model.sigma);
    K = model.kernelFunc(X,model.Xtrain,model.kernelArgs);
    N = (1/model.nTrain)*ones(model.nTrain);
    K = K - N*K - K*N + N*K*N;
    Xreduced = K*model.basis; % Note: basis is orthogonal so inverse is transpose
end