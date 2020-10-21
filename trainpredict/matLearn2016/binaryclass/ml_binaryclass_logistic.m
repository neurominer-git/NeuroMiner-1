function [model] = ml_binaryclass_logistic(X,y,options)
% ml_binaryclass_logistic(X,y,options)
%
% Description:
%	 - Fits a classifier using logistic regression
%
% Options:
%    - addBias: adds a bias variable (default: 1)
%    - lambdaL2: strenght of L2-regularization parameter (default: 0)
%
% Authors:
% 	 - Neil Newman (2014)

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[addBias,lambdaL2] = myProcessOptions(options,'addBias',1,'lambdaL2',0);

% Add optional bias
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Try to optimize
optimOptions.Display = 0;
optimOptions.useMex = 1;

if lambdaL2 == 0
    w = minFunc(@logisticLoss,randn(nFeatures,1),optimOptions,X,y);
else
    w = minFunc(@logisticLossL2,randn(nFeatures,1),optimOptions,X,y,lambdaL2);
end

% Model outputs
model.name = 'Logistic Regression';
model.addBias = addBias;
model.w = w;
model.predict = @predict;
model.predictProb = @predictProb;
end

function p = predict(model,Xhat)
% Prediction funcion
[nTest,nFeatures] = size(Xhat);
p.yhat = zeros(nTest, 1);
probs = predictProb(model, Xhat);
p.D = probs(:, 2)-0.5;

% Threshold at 0.5
index = p.D >= 0;
p.yhat = sign(index - 0.5);
end

function [probs] = predictProb(model,Xhat)
% Returns a vector of entries [p(y=-1), p(y=1)]
[nTest,nFeatures] = size(Xhat);
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end
prob = sigmoid(Xhat * model.w);
probs = [1 - prob, prob];
end

function [f,g] = logisticLoss(w,X,y)
% Logistic loss function with training weights
[nTrain,nFeatures] = size(X);
yXw = y.*(X*w);
f = sum(log(1 + exp(-yXw)));
if isinf(f)
    % Fallback on log-sum-exp trick if you overflow
    f = sum(logsumexp([zeros(nTrain,1) -yXw]));
end
g = -(X.'*(y./(1+exp(yXw))));
end

function [f,g] = logisticLossL2(w,X,y,lambda)
% Logistic loss function with lambda and training weights
yXw = y.*(X*w);
f = sum(log(1 + exp(-yXw)));
if isinf(f)
    % Fallback on log-sum-exp trick if you overflow
    f = sum(logsumexp([zeros(n,1) -yXw]));
end
f = f + (lambda/2)*(w'*w);
g = -(X.'*(y./(1+exp(yXw)))) + lambda*w;
end

function [lse] = logsumexp(b)
% Compute log(sum(exp)) across columns without overflowing
B = max(b,[],2);
lse = log(sum(exp(b-repmat(B,[1 size(b,2)])),2))+B;
end

function [val] = sigmoid(z)
% Sigmoid function
val = 1 ./ (1 + exp(-z));
end