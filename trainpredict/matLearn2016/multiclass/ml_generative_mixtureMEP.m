function [model] = ml_generative_mixtureMEP(X, ~, options)
% ml_generative_mixtureMEP(X, ~, options)
%
% Description:
%	 - Finds the best fit for a homogenous (fixed exponent) Multivariate 
%      Exponential Power mixture model. Set kappa = 1 for multivariate 
%      Laplacian or kappa = 2 for multivariate Gaussian.
%
% Options:
%   - init: Character string. 'rand' for random initialization, 'kpp' for 
%           k-means++ initialization
%   - nMixtures: Integer representing the target number of clusters
%   - maxIter: Number of times one run of the EM algorithm should iterate
%   - optTol: tuning parameter for EM algorithm defining the minimum change
%             from epoch to epoch
%   - nRestarts: how many times EM algorithm runs
%
% Author:
% 	 - Geoffrey Roeder (2016), based on MSc thesis of Sascha Brauer (2014)

[nTrain,nFeatures] = size(X);
[nMixtures, init, kappa, maxIter, optTol, nRestarts] = ...
myProcessOptions(options, 'nMixtures', 3, 'init', 'rand', ...
                'kappa', 1, ... % Multivariate Laplacian mixture
                'maxIter', 1000, ...
                'optTol', 1e-14, ...
                'nRestarts', 5);

likelihood = -inf;
model = [];
Xi = zeros(nTrain,nMixtures);
pnc = zeros(nTrain,nMixtures);

% ensure iid assumption holds
X = X(randperm(nTrain), :);

for i = 1:nRestarts+1
    % Initialize by giving each cluster a random assignment of points.
    [mu, sigma, w] = initialize_params(X, nMixtures, init);

    % Iterate EM algorithm until convergence
    previousLikelihood = -inf;
    for j = 1:maxIter
        % E-step-----------------------------------------------------------
        % find the Xi matrix, exponents of the MEP pdf
        for c = 1:nMixtures
           Xc = bsxfun(@minus, X, mu(c,:));
           % Fast method using Cholesky
           [R, err] = cholcov(sigma(:,:,c),0);
           if err ~= 0
               error('Covariance matrix was positive semi-definite')
           end
           Xi(:,c) = sum((Xc / R).^2,2);
        end
        
        % evaluate responsibilities for each point
        for c = 1:nMixtures 
            pnc(:,c) = simple_meppdf(X, mu(c,:), sigma(:,:,c), ...
                                     kappa) * w(c); 
        end
        sum_p = sum(pnc, 2);
        rnc = bsxfun(@rdivide,pnc,sum_p);
        rnc(isnan(rnc)) = 0;
        currentLikelihood = sum(log(sum_p));
        
        % M-step ----------------------------------------------------------
        % update mixture weights
        w = mean(rnc);
       for c = 1:nMixtures
           Xc = bsxfun(@minus, X, mu(c,:));
           U = Xi(:,c).^(kappa/2-1);
           U(isinf(U)) = 0;
           U(isnan(U)) = 0;
           rXi_c = (rnc(:,c).*U)';
           mu(c,:) = rXi_c*X / sum(rXi_c);
           rXiXc = bsxfun(@times, rnc(:,c), bsxfun(@times, ...
                                             U, Xc));
           sig_numer = rXiXc'*Xc;                       
           sig_denom = sum(rnc(:,c));
           sigma(:,:,c) = (kappa/2).*sig_numer./sig_denom + ...
                          1e-8*eye(nFeatures);
       end
        currentModel.w = w;
        currentModel.mu = mu;
        currentModel.sigma = sigma;
        
        % Check optimality condition --------------------------------------
        if currentLikelihood <= previousLikelihood + optTol
            break;
        else
            previousLikelihood = currentLikelihood;
        end
    end
    
    if currentLikelihood > likelihood
        likelihood = currentLikelihood;
        model = currentModel;
    end
end

if likelihood == -inf
    error('Failed to fit MEP Mixture Model. Try different parameters.')
end

model.name = 'Generative MEP Mixture';
model.nMixtures = nMixtures;
model.kappa = kappa;
model.predict = @predict;
end

function [y] = simple_meppdf(X, mu, Sigma, kappa)
[~ ,d] = size(X);

% Center data
X_cent = bsxfun(@minus, X, mu);

% Evaluate Cholesky decomp
[R, e] = cholcov(Sigma,0);
% Check for positive definite versus positive semi-definite
if e > 0
    error('Bad Sigma: was not positive definite');
end
% Standardize centered data
Xs = X_cent / R;
log_sqr_det_sig = sum(log(diag(R)));
exponent = sum(Xs.^2, 2).^(kappa / 2);
C_exp = 1 + d / kappa;
C = d * gamma(d / 2) ./ ((pi^(d/2))*gamma(C_exp)*2^(C_exp));
y = exp( -0.5 * exponent - log_sqr_det_sig + log(C));    
end

function [mu, sigma, w] = initialize_params(X, nMixtures, init)
    [nTrain, nFeatures] = size(X);   
    if strcmp(init, 'kpp')
        options = [];
        options.k = nMixtures;
        options.kpp = 1;
        options.plot = 0;
        model = ml_unsupervised_clusterKmeans(X, options);
        inds = model.assignment;
    elseif strcmp(init, 'rand')
       % generate random assignments for points
       % rng(1)
       inds = mod(randperm(nTrain), nMixtures) + 1;
    else
       error('Unrecognized initialization type: %s', init);
    end
        
    mu = zeros(nMixtures,nFeatures);
    sigma = zeros(nFeatures,nFeatures,nMixtures);
    w = zeros(nMixtures);
    for c = 1:nMixtures
        X_c = X(inds==c, :);
        % Fit a Gaussian to the cluster
        mu(c,:) = mean(X_c,1);
        sigma(:,:,c) = cov(X_c) + 1e-8*eye(nFeatures); % make pos. def.
        % initialize mixture weights 
        w(c) = sum(inds == c) / nTrain;
    end
end

function [density] = predict(model,Xhat)
    [nTest, ~] = size(Xhat);
    w = model.w;
    mu = model.mu;
    sigma = model.sigma;
    kappa = model.kappa;
    nMixtures = model.nMixtures;
    density = zeros(nTest, 1);
    for c = 1:nMixtures
        density = density + ...
                  simple_meppdf(Xhat, mu(c,:), sigma(:,:,c),kappa) * w(c);
    end
end
