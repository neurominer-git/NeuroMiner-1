function [model] = ml_regression_leastSquares(X,y,options)
% w = ml_regression_leastSquares(X,y)
%
% Returns Least Squares Solution:
%  min_w sum (Xw - y).^2 by solving the Normal Equations
%
% Options:
%   lambdaL2: scalar regularizer value for loss function
%   weights: n x n matrix of weights
%
% Authors:
%   - Mark Schmidt (2014); adapted for matLearn by Geoffrey Roeder (2016)

if nargin < 3
    options = [];
end

model.name = 'Least Squares';

[lambdaL2,weights] = myProcessOptions(options,'lambdaL2',0,'weights',[]);
[nInstances,nVars] = size(X);

if isempty(weights)
    if lambdaL2 == 0
        % Least Squares
        w = (X'*X)\(X'*y);
    else
        model.name = 'Ridge Regression';
        % Least Squares w/ L2 Prior (Ridge Regression)
        w = (X'*X + lambdaL2*eye(nVars))\(X'*y);
    end
else
    W = diag(weights);
    if lambdaL2 == 0
        model.name = 'Weighted Least Squares';
        w = (X'*W*X)\(X'*W*y);
    else
        model.name = 'Weighted Ridge Regression';
        w = (X'*W*X + lambdaL2*eye(nVars))\(X'*W*y); 
    end
end

model.weights = w;
model.predict = @(model,X)X*model.weights;
model.errL2Func = @(model,X,y)sum((X*model.weights-y).^2);
model.lossFunc = @(model,X,y)sum((X*model.weights-y).^2);
