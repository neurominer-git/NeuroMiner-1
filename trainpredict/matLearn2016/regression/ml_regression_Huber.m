function [ model ] = ml_regression_Huber(X,y,options)
% ml_regression_Huber(X,y,options)
%
% Description:
%	 - Predicts parameters based on Huber loss function
%       (combines both l1 and l2 losses)
%
% Options:
%    - weights: weight on each data point(default: 0)
%    - lambdaL2: strength of L2-regularization parameter (default:0)
%    - addBias: adds a bias variable (default: 0)
%    - epsilon: boundary between minimizing either
%       Least Sqaures (LS) a Least Absolute Deviations (LAD)
%       residual values < transition uses LS (default: 0.9)
%
% Author: Adrian Wong (2014)

[nTrain,nFeatures] = size(X);

% Use default options if not user specified
[z,lambdaL2,addBias,epsilon] = ...
    myProcessOptions(options,'weights',ones(nTrain, 1),'lambdaL2',0, ...
                     'addBias',1,'epsilon',0.9);

% Add bias if required
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Optimization options
optimOptions.Display = 0;
optimOptions.useMex = 1;

if lambdaL2 == 0
    % Find regression weights (without L2 regularization)
    w = minFunc(@HuberLoss,randn(nFeatures,1),optimOptions,X,y,epsilon,z);
else
    % Find regression weights with L2 regularization
    w = minFunc(@HuberLossL2,randn(nFeatures,1),optimOptions,X,y,epsilon,z,lambdaL2);
end

% Model outputs
model.name = ['Huber Loss Linear Regression \epsilon=',num2str(epsilon)];
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

function [f,g] = HuberLoss(w,X,y,epsilon,z)
% Huber loss function without L2 regularization
r = X*w-y;
absr = abs(r);
id = abs(r) <= epsilon;
id_2 = r < 0;
f = z'*(id.*(0.5*absr.^2)+(~id).*(epsilon*absr-0.5*epsilon^2));
g = X'*diag(z)*diag(id)*(X*w-y) + (-epsilon)*X'*diag(z)*id_2+ epsilon*X'*diag(z)*(~id_2);
end

function [f,g] = HuberLossL2(w,X,y,t,z,lambda)
% Huber loss function with L2 regularization
r = X*w-y;
absr = abs(r);
id = absr <= t;
id_2 = r < 0;
f = z'*(id.*(0.5*absr.^2)+(~id).*(t*absr-0.5*t^2)) + lambda*(w'*w);
g = X'*diag(z)*diag(id)*(X*w-y) + (-t)*X'*diag(z)*diag(id_2)+ t*X'*diag(z)*diag(~id_2) + 2*lambda*w;
end