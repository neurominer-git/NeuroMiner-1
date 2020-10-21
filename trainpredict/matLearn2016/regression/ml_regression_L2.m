function [model] = ml_regression_L2(X,y,options)
% ml_regression_L2(X,y,options)
%
% Description:
%	 - Fits a linear regression model by minimizing the squared loss
%
% Options:
%	 - weights: vector of training data weights (default: vector of 1's)
%	 - lambdaL2: strength of L2-regularization parameter (default: 0)
%	 - addBias: accepts 0 or 1. If 1, adds bias to X (default: 1)
%	 - features: indices for one-of-k regression onto subset of features
%
% Authors:
% 	- Yan Zhao (2014); one-of-k regresssion added by Geoffrey Roeder (2016)

[nTrain,nFeatures] = size(X);

% Generate default model options if not given by input
[z,lambdaL2,addBias,features] = myProcessOptions(options,'weights', ...
                                        ones(nTrain,1),'lambdaL2',0, ...
                                        'addBias',1,'features',1:nFeatures);
X = X(:,features);
[nTrain,nFeatures] = size(X);

% Add bias variable if required
if addBias
    X = [ones(nTrain,1) X];
    features = [1 features+1];
    nFeatures = nFeatures + 1;
end

% Optimization options for minFunc
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Compute regression weights using minFunc
w = minFunc(@squaredLossL2,randn(nFeatures,1),optimOptions,X,y,lambdaL2,z);

% Model outputs
model.name = 'Squared Loss Linear Regression';
model.w = w;
model.addBias = addBias;
model.predict = @predict;
model.features = features;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
% Add optional bias variable
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
    nFeatures = nFeatures + 1;
end
% select one-of-k subset
Xhat = Xhat(:,model.features);
yhat = Xhat*model.w;
end

function [f,g] = squaredLossL2(w,X,y,lambda,z)
% Squared loss function with training sample weights and L2 regularization
f = (X*w-y)'*diag(z)*(X*w-y) + lambda*(w'*w);
g = 2*X'*diag(z)*(X*w-y) + 2*lambda*w;
end
