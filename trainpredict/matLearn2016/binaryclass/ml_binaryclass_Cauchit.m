function [ model ] = ml_binaryclass_Cauchit(X,y,options)
% ml_binaryclass_Cauchit(X,y,options)
%
% Description:
%	 - Classifies data using generalized linear model with a Cauchit link
%	 - In order to fit the data, link is passed to
%       MATLAB's glmfit, along with the design matrix and label vector
%
% Options:
%    - dist: distribution (default: binomial)
%	 - thresh: classification threshold (default: 0.5)
%
% Authors:
% 	 - Delaram Behnami (2014)

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[dist, thresh] = myProcessOptions(options,'dist','binomial','thresh',0.5);

% Construct Cauchit link by entering definition, derivative and inverse
Cauchit.Link = @(X) tan(pi*(X-.5));
Cauchit.Derivative = @(X) pi*(csc(X.*pi)).^2;
Cauchit.Inverse = @(X) (2*atan(X)+pi)/(2*pi);

% Fit generalized linear model to data using extreme link
% Obtain fitting vector of weights (w)
y = (y + 1)./2;
w = glmfit(X, [y ones(nTrain,1)], dist, 'link', Cauchit);

% Construct fitting model
model.w = w;
model.link = Cauchit;
model.predict = @predict;
model.predictProb = @predictProb;
model.name = 'Cauchit Loss Binary Classification';
model.thresh = thresh;
end

function p = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Obtain predicted regression probability
probs = predictProb(model,Xhat);
p.prob = probs(:,2);

% Predict class by thresholding the probability
index = p.prob >= model.thresh;
p.yhat = sign(index - 0.5);
end

function [probs] = predictProb(model, Xhat)
% Returns a vector of entries [p(y=-1), p(y=1)]
prob = glmval(model.w,Xhat,model.link);
probs = [1-prob, prob];
end