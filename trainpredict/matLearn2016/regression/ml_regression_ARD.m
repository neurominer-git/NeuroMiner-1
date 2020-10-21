function [model] = ml_regression_ARD(X, y, options)
% ml_regression_ARD(X, y, options)
%
% Description:
%	 - Performs type 2 maximum likelihood for linear regression using
%          Mackay's implementation
%
% Options:
%    - sigma2: variance parameter (default: 1e-5)
%    - addBias: whether or not to add an intercept (default: 1)
%    - scale: whether or not to standardize columns of feature variables
%       (default: 0)
%
% Authors:
%    - Ken Lau (2014)

[nTrain,nFeatures] = size(X);

[sigma2, addBias, scale] = myProcessOptions(options, 'sigma2', 1e-5,...
                                                     'addBias', 1,...
                                                     'scale', 0);

% Scale the feature variable matrix X if required
if scale
    [X, muX, sX] = standardizeCols(X);
    model.muX = muX;
    model.sX = sX;
end

% Add bias variable if required
if addBias
    X = [ones(nTrain,1), X];
    nFeatures = nFeatures + 1;
end

% Initalize hyper-parameter to optimize
lambda = ones(nFeatures, 1);

% Optimize hyper-parameter lambda
while 1
    lambdaOld = lambda;
    
    C = diag(sigma2*ones(nTrain,1)) + X*diag(lambda)*X';
    w = diag(lambda)*X'*(C\y);
    V = diag(lambda) - diag(lambda)*X'*inv(C)*X*diag(lambda);
    
    lambda = w.^2./(1-diag(V)./lambda);
    lambda(isinf(lambda)) = 1e-15;
    
    if max(abs(lambdaOld - lambda)) < 1e-3
        break;
    end
end

% Returns the index of features with non-zero weights
idxFeat = 1:nFeatures;
% Obtain indices of non-zero weights
expr = abs(w) > 1e-12;
if addBias
    % Get rid of the intercept term
    idxFeat = idxFeat(1:end-1);
    expr = expr(2:end);
    idxFeat = idxFeat(expr);
else
    idxFeat = idxFeat(expr);
end

% Model outputs
model.name = 'ARD';
model.addBias = addBias;
model.scale = scale;
model.w = w;
model.idxFeat = idxFeat;
model.predict = @predict;
model.featuresSelected = idxFeat;

end

function [yhat] = predict(model, Xhat)
% Prediction function
% Back transformation if scaled
if model.scale
    Xhat = standardizeCols(Xhat, model.muX, model.sX);
end

% Add on the bias term if biasis included
if model.addBias
    Xhat = [ones(size(Xhat,1),1), Xhat];
end
yhat = Xhat*model.w;
end
