function model = ml_binaryclass_elastic_net(X,y,options)
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
% 	 - Nikolaos Koutsouleris (2017)

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[addBias,lambdaL1, lambdaL2] = myProcessOptions(options,'addBias',1,'lambdaL1',0, 'lambdaL2',0);

% Add optional bias and prepare for optimization
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1; 
end

lambdaL1 = lambdaL1*ones(nFeatures,1);
w_init = zeros(nFeatures,1);

if addBias, lambdaL1(1)=0; end

% Optimize
funObj = @(w)LogisticLoss(w,X,y);
funObjL2 = @(w)penalizedL2(w,funObj,lambdaL2);
w = L1General2_PSSgb(funObjL2,w_init,lambdaL1);

% Model outputs
model.name = 'Elastic Net Logistic Regression';
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
try
    prob = sigmoid(Xhat * model.w);
catch
    fprintf('problem')
end
probs = [1 - prob, prob];
end

function [val] = sigmoid(z)
% Sigmoid function
val = 1 ./ (1 + exp(-z));
end