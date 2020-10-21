function [model] = ml_regression_NB(X,y,options)
% ml_regression_NB(X,y,options)
%
% Description:
%	 - Fits a linear regression model by minimizing the "Naive Bayes" 
%      squared loss, where each feature is assumed to be independent of all
%      the others
%
% Options:
%	 - weights: vector of training data weights (default: vector of 1's)
%	 - lambdaL2: strength of L2-regularization parameter (default: 0)
%	 - addBias: accepts 0 or 1. If 1, adds bias to X (default: 1)
%
% Authors:
%    - Scott Sallinen (2014)

[nTrain,nFeatures] = size(X);

% Generate default model options if not given by input
[z,lambdaL2,addBias] = myProcessOptions(options, 'weights', ...
                                        ones(nTrain,1), 'lambdaL2', 0, ...
                                        'addBias', 1);

% Add a bias column to X if required
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Optimization options
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Compute regression weights using minFunc
w = minFunc(@squaredLossNB, randn(nFeatures, 1), optimOptions, X, y, ...
            lambdaL2, z);

% Model outputs
model.name = 'NB Squared Loss Linear Regression';
model.w = w;
model.addBias = addBias;
model.predict = @predict;
end


function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
% Add bias if required
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end
yhat = Xhat*model.w;
end

function [f,g] = squaredLossNB(w,X,y,lambda,z)
% Naive Bayes squared loss function
[nTrain, nFeatures] = size(X);
f = lambda*(w'*w);
g = 2*lambda*w;
for j=1:nFeatures
    % Naive Bayes summation
    difj = X(:,j)*w(j)-y;
    f = f + difj'*diag(z)*difj;
    g(j) = 2*X(:,j)'*diag(z)*difj;
end
end