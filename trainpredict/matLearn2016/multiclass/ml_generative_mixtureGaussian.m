function [model] = ml_generative_mixtureGaussian(X, ~, options)
% ml_multiclass_generativeMixtureGaussian(X,~,options)
%
% Description:
%	 - Finds the best fit of k Gaussians to the data using Expectation
%       Maximization. Targets y are expected to be empty and are ignored. The
%       prediction function returns an [nTest, 1] vector with the likelihood of
%       each data point in the test set.
%
% Options:
% 	 - nMixtures: the number of Gaussians to fit to the data (default: 2)
%	 - maxIter: maximum number of iterations before terminating (default: 100)
%	 - optTol: tolerance of data log likelihood to accept convergence
%       	   (default: 1e-14)
%	 - nRestarts: number of random restarts allowed. (default: 5)
%
% Authors:
% 	 - Neil Traft (2014)

[nTrain,nFeatures] = size(X);
[nMixtures, maxIter, optTol, nRestarts] = myProcessOptions(options,...
    'nMixtures',2,...
    'maxIter',1000,...
    'optTol',1e-14,...
    'nRestarts',5);

lik = -inf;
model = [];
for i = 1:nRestarts+1
    % INITIALIZATION PHASE -------------------
    w = zeros(nTrain, nMixtures);
    p = ones(nMixtures, 1)/nMixtures;
    
    % Initialize by giving each cluster a random assignment of points.
    ind = randperm(nTrain);
    n = nTrain/nMixtures;
    X = X(randperm(nTrain), :);
    % generate random assignments for points
    inds = mod(randperm(nTrain), nMixtures) + 1;
    
    for c = 1:nMixtures
        % Assign points
        rows = floor((c-1)*n+1) : floor(c*n);
        %subX = X(ind(rows), :);
        subX = X(inds==c, :);
        % Fit a Gaussian to the cluster
        mu(:,c) = mean(subX,1);
        sigma(:,:,c) = cov(subX) + 1e-8*eye(nFeatures); % make pos. def.
    end
    
    % E-STEP PHASE ------------------------
    % Iterate until convergence.
    prevLik = -inf;
    for j = 1:maxIter
        % E step: Perform soft cluster assignment
        for c = 1:nMixtures
            w(:,c) = mvnpdf(X, mu(:,c)', sigma(:,:,c)) * p(c);
        end
        
        summ = sum(w,2);
        w = bsxfun(@rdivide, w, summ);
        tempLik = sum(log(summ));
        
        % M step: Optimize Q
        p = mean(w);
        sqrtR = sqrt(w);
        
       for c = 1:nMixtures
            Xnorm = bsxfun(@minus, X, mu(:,c)');
            mu(:,c) = (w(:,c)'*X) / sum(w(:,c));
            rX = bsxfun(@times, sqrtR(:,c), Xnorm);
            scatter = rX'*rX;
            sigma(:,:,c) = scatter / sum(w(:,c)) + 1e-6*eye(nFeatures);
       end
        
        % Save and (possibly) visualize.
        tempModel.p = p;
        tempModel.mu = mu;
        tempModel.sigma = sigma;
        
        % Test if we've reached an optimum (the objective is no longer moving).
        if tempLik <= prevLik + optTol
            break;
        else
            prevLik = tempLik;
        end
        
    end
    
    if tempLik > lik
        lik = tempLik;
        model = tempModel;
    end
end

model.name = 'Generative Gaussian Mixture Model';
model.nMixtures = nMixtures;
model.predict = @predict;
end

function [lik] = predict(model,Xhat)
[nTest, nFeatures] = size(Xhat);
p = model.p;
mu = model.mu;
sigma = model.sigma;
nMixtures = model.nMixtures;
lik = zeros(nTest, 1);
for c = 1:nMixtures
    % r_ic = p(x_i | \mu_c, \sigma_c) \pi_c
    lik = lik + mvnpdf(Xhat, mu(:,c)', sigma(:,:,c)) * p(c);
end
end
