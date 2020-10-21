function [ model ] = ml_binaryclass_extreme(X,y,options)
% ml_binaryclass_extreme(X,y,options)
%
% Description:
%	 - Classifies data using generalized linear model with a extreme link
%	 - In order to fit the data, link is passed to
%       MATLAB's glmfit, along with the design matrix and label vector
%
% Options:
%    - dist: distribution (default: binomial)
%    - thresh: classification threshold (default: 0.5)
%
% Authors:
% 	 - Shekoofeh Azizi (2014)

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[dist, thresh] = myProcessOptions(options,'dist','binomial','thresh',0.5);

% Construct extreme link by entering definition, derivative and inverse
extreme.Link = @(x) -log(-log(x));
extreme.Derivative =@(x) -1./(x.*log(x));
extreme.Inverse = @(x) exp(-exp(-x));

% Fit generalized linear model to data using extreme link
% Obtain fitting vector of weights (w)
y = (y + 1)./2;
w = glmfit(X, [y ones(nTrain,1)], dist, 'link', extreme);

% Construct fitting model
model.w = w;
model.link = extreme;
model.predict = @predict;
model.predictProb = @predictProb;
model.name = 'Extreme Loss Binary Classification';
model.thresh = thresh;
end

function p = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Obtain predicted regression probability
p.yhat = zeros(nTest, 1);
probs = predictProb(model, Xhat);
p.D = probs(:, 2)-0.5;

% Threshold at 0.5
index = p.D >= 0;
p.yhat = sign(index - 0.5);
end

function [probs] = predictProb(model, Xhat)
% Returns a vector of entries [p(y=-1), p(y=1)]
prob = glmval(model.w,Xhat,model.link);
probs = [1-prob, prob];
end