function [model] = ml_binaryclass_probit(X,y,options)
% ml_binaryclass_probit(X,y,options)
%
% Description:
%	 - Fits a probit regression classifier by minimizing the negative likelihood
%
% Options:
%    - addBias: adds a bias variable (default: 1)
%    - lambdaL1: strength of L1-regularization parameter (default: 0)
%    - lambdaL2: strength of L2-regularization parameter (default: 0)
%
% Authors:
% 	 - Rindra Ramamonjison (2014)
%    - Mark Schmidt
%
% Source:
%    - Mark Schmidt, Glenn Fung, Romer Rosales. Optimization Methods for
%       L1-Regularization. UBC Technical Report TR-2009-19, 2009

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[addBias,lambdaL1, lambdaL2] = myProcessOptions(options,'addBias',1,'lambdaL1',0,'lambdaL2',0);

% Add optional bias
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Try to optimize
optimOptions.Display = 0;
optimOptions.useMex = 1;

if lambdaL2 == 0 && lambdaL1 == 0
    % No regularization
    w = minFunc(@probitLoss,randn(nFeatures,1),optimOptions,X,y);
elseif lambdaL1 == 0
    % L2 regularization
    w = minFunc(@probitLossL2,randn(nFeatures,1),optimOptions,X,y,lambdaL2);
elseif lambdaL2 == 0
    % L1 regularization
    % Set penalty weights
    lambda = lambdaL1*ones(nFeatures,1);
    % Objective function handle
    funObj = @(w)probitLoss(w,X,y);
    % L1 augmented objective function
    w = L1General2_PSSgb(funObj,zeros(nFeatures,1),lambda);
else
    % Elastic net regularization
    % Set penalty weights
    lambda = lambdaL1*ones(nFeatures,1);
    % Objective function handle
    funObj = @(w)probitLossL2(w,X,y, lambdaL2);
    % Regularized objective function
    w = L1General2_PSSgb(funObj,zeros(nFeatures,1),lambda);
end

% Model outputs
model.name = 'Probit Loss Binary Classification';
model.addBias = addBias;
model.w = w;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end

prob = normcdf(Xhat*model.w);

% Threshold at 0.5
index = prob >= 0.5;
yhat = sign(index - 0.5);
end

function [f,g] = probitLoss(w,X,y)
% Probit loss function with training weights
Xw = X*w;
% Transform 0/1 encoding of y to -1/1
s = y;
if sum(y == -1) == 0
    s = 2*y-1;
end
probit = normcdf(s.*Xw);
f = - sum(log(probit));

c = s.*normpdf(s.*(Xw))./probit;
g = - X'*c;
end

function [f,g] = probitLossL2(w,X,y,lambda)
% Probit loss with training weights and lambda
Xw = X*w;
% Transform 0/1 encoding of y to -1/1
s = y;
if sum(y == -1) == 0
    s = 2*y-1;
end
probit = normcdf(s .* Xw);
f = - sum(log(probit)) + (lambda/2)*(w'*w);

c = s.*normpdf(s.*(Xw))./probit;
g = - X'*c + lambda*w;
end