function [model] = ml_ordinal_logistic(X,y,options)
% ml_ordinal_logistic(X,y,options)
%
% Description:
%	 - Classification using ordinal logistic regression
%
% Options:
%    - verbose: verbosity of optimization subroutine (default: 1)
%    - lambdaL2: strength of L2-regularization parameter (default: 0)
%
% Authors:
% 	 - Mark Schmidt (2014)


[n,p] = size(X);

if nargin < 3
    options = [];
end

[verbose,lambdaL2] = myProcessOptions(options,'verbose',1,'lambdaL2',0);
k = options.nClasses;

if verbose
    optimoptions.verbose = 3;
else
    optimoptions.verbose = 0;
end

w = zeros(p,1);
gamma = ones(k-2,1);
LB = [-inf(p,1);zeros(k-2,1)];
UB = inf(p+k-2,1);
if lambdaL2 == 0
    model.name = 'Ordinal Logistic Regression';
    funObj = @(w)ml_ordinal_logistic_loss(w,X,y,k);
else
    model.name = 'Ordinal Logistic Regression with L2-Regularization';
    funObj_sub = @(w)ml_ordinal_logistic_loss(w,X,y,k);
    % First element is assumed to be the bias
    lambdaVect = [0;lambdaL2*ones(p-1,1);zeros(k-2,1)]; 
    funObj = @(w)ml_penalized_L2(w,funObj_sub,lambdaVect);
end

w_gamma = minConf_TMP(funObj,[w(:);gamma(:)],LB,UB,optimoptions);
w = w_gamma(1:p);
gamma = [-inf;0;cumsum(w_gamma(p+1:end));inf];

model.nClasses = k;
model.weights = w;
model.thresh = gamma;
model.predict = @predict;
model.lossFunc = @loss;

end

function y = predict(model,X)
k = model.nClasses;
w = model.weights;
gamma = model.thresh;
[n,p] = size(X);
z = X*w;
y = zeros(n,1);
for c = 1:k
    y(z > gamma(c)) = c;
end
end

function f = loss(model,X,y)
w = model.weights;
gamma = model.thresh;
sigmoid1 = 1./(1+exp(X*w - gamma(y+1)));
sigmoid2 = 1./(1+exp(X*w - gamma(y)));
f = -sum(log(sigmoid1-sigmoid2+eps))';
end
