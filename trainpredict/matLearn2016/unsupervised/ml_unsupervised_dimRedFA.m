function [model] = ml_unsupervised_dimRedFA(X, options)
% ml_unsupervised_dimRedFA(X,options)
%
% Description:
%	 - Fit factor analysis model to dataset assuming nComponents factors in
%	   order to reduce dimensionality of input data X
%
% Options:
%   - nComponents: target number of clusters 
%   - doPPCA: assume all factors have a shared, single variance s^2 so
%             that the covariance matrix is s^2*eye(nVars)
%
% Authors:
% 	 - Mark Schmidt (2014)

if nargin < 2
    options = [];
end

[nInstances,nVars] = size(X);

[nComp, doPPCA] = myProcessOptions(options, 'nComponents', ...
                                   nVars, 'doPPCA', 0);

% Standardize Columns
[X,mu,sigma_data] = standardizeCols(X);

% Form empirical covariance
C = (1/nInstances)*X'*X;

% First fit regular PCA to initialize basis
[U,S,V] = svd(C);
W = V(:,1:nComp);

% Now fit probabilistic PCA to initialize variance
sigma = 1;
funObj = @(Wsigma)ml_ppca_loss(Wsigma,nVars,nComp,C);
options.numDiff = 1;
options.maxIter = 10000;
options.optTol = 1e-16;
options.verbose = 0;
LB = [-inf(nVars*nComp,1);0];
UB = inf(nVars*nComp+1,1);
Wsigma = minConf_TMP(funObj,[W(:);sigma],LB,UB,options);
W = reshape(Wsigma(1:nVars*nComp),nVars,nComp);
sigma = Wsigma(end);

% Now fit factor analysis
if ~doPPCA
    sigma = sigma*ones(nVars,1);
    funObj = @(Wsigma)ml_factor_analysis_loss(Wsigma,nVars,nComp,C);
    LB = [-inf(nVars*nComp,1);zeros(nVars,1)];
    UB = inf(nVars*nComp+nVars,1);
    Wsigma = minConf_TMP(funObj,[W(:);sigma],LB,UB,options);
    W = reshape(Wsigma(1:nVars*nComp),nVars,nComp);
    sigma = Wsigma(nVars*nComp+1:end);
end 

model.mu = mu;
model.W = W;
model.sigma = sigma;
model.sigma_data = sigma_data;
model.nComponents = nComp;
model.reduceFunc = @reduce;
end

function [Xreduced] = reduce(model,X)
    X = standardizeCols(X,model.mu,model.sigma_data);
    Xreduced = (model.W\X')'; % Note: basis may not be orthogonal, use pinv
end
