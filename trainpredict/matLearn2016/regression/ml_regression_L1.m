function [model] = ml_regression_L1(X,y,options)
% ml_regression_L1(X,y,options)
%
% Description:
%    - Learns a set of weights for linear regression which minimizesthe
%       absolute (L1) loss
%
% Options:
%  	 - weights: vector of training data weights (default: vector of 1's)
%	 - lambdaL2: strength of L2-regularization parameter (default: 0)
%	 - addBias: accepts 0 or 1. If 1, adds bias to X (default: 1)
%
% Authors:

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[z,lambdaL2,addBias] = myProcessOptions(options,'weights',ones(nTrain,1),'lambdaL2',0, 'addBias', 1);

% Add bias variable if necessary
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Inputs to linprog/quadprog
% Minimize \sigma_{i=1}^N v_i (+ lambdaL2 * L2 regularization)
%   subject to constraints: X_i'w - v_i <= y_i and -X_i'w - v_i <= -y_i
f = [zeros(nFeatures,1); z];
b = [y;-y];
A = [X -eye(nTrain);-X -eye(nTrain)];
options_ = optimset('Display', 'off');
% Linear programming without L2 regularization
if lambdaL2 == 0
    [wv] = linprog(f,A,b,[],[],[],[],[], options_);
% Quadratic programming with L2 regularization
else
    H = 2*lambdaL2*[eye(nFeatures) zeros(nFeatures, nTrain); zeros(nTrain, nFeatures), zeros(nTrain, nTrain)];
    [wv] = quadprog(H,f,A,b,[],[],[],[],[], options_);
end
% Model ouputs
model.w = wv(1:nFeatures);
model.name = 'Absolute Loss Linear Regression';
model.addBias = addBias;
model.predict = @predict;

end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
% Add optional bias variable
if model.addBias
    Xhat = [ones(nTest, 1), Xhat];
    nFeatures = nFeatures + 1;
end
yhat = Xhat*model.w;
end