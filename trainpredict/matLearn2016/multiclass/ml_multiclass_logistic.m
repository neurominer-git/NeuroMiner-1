function [model] = ml_multiclass_logistic(X, y, options)
% ml_multiclass_logistic(X,y,options)
%
% Description:
%	 - Classification using multinomial logistic regression
%
% Options:
%    - addBias: adds a bias variable (default: 1)
%    - standardize: standardize the data (default: 0)
%    - lambdaL2: strength of L2-regularization parameter (default: 0)
%
% Authors:
% 	 - Keyulu Xu (2014)

[nTrain,nFeatures] = size(X);

% Default outputs
[addBias,standardize,lambdaL2] = myProcessOptions(options,'addBias',1,...
                                                  'standardize',0,...
                                                  'lambdaL2',0);

% Calculate number of unique classes in y
classes = unique(y);
nClasses = length(classes);

% Transform the input data X
if standardize
    mu = mean(X);
    sigma = std(X);
    X = X - repmat(mu, nTrain, 1);
    X = X ./ repmat(sigma, nTrain, 1);
    
    % Make sure we can apply same transformation at testing time
    model.standardize = 1;
    model.mu = mu;
    model.sigma = sigma;
end

% Optionally add bias
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Optimization parameters
optimOptions.Display = 0;
optimOptions.useMex = 0;

if lambdaL2 == 0
    % Without regularization
    w = minFunc(@logisticLoss,randn(nFeatures*nClasses,1),optimOptions,...
                X,y,classes);
else
    w = minFunc(@logisticLossL2,randn(nFeatures*nClasses,1),...
                optimOptions,X,y,lambdaL2,classes);
end

w = reshape(w, nFeatures, nClasses);

% Model outputs
model.w = w;
model.addBias = addBias;
model.standardize = standardize;
model.name = 'Multiclass Logistic Classification';
model.nClasses = nClasses;
model.classes = classes;
model.predict = @predict;
end

function [f, g] = logisticLoss(w, X, y, classes)
% Logistic loss function
% w is nFeatures * nClasses
% X is nTrain * nFeatures
% y is nTrain * nClasses
% z is nTrain * 1
[~, nFeatures] = size(X);
nClasses = length(classes);

% Reshape so that we can multiply
w = reshape(w, [nFeatures, nClasses]); 

yBin = bsxfun(@eq, y, classes');
tmp = X*w - repmat(logsumexp(X*w), 1, nClasses);
f = -sum(sum(tmp .* yBin, 2));

% Compute gradient
mu = exp(tmp);
g = -(X'*(yBin - mu));

% Reshape for optimization
g = reshape(g,[nFeatures * nClasses, 1]); 
end

function [f, g] = logisticLossL2(w,X,y,lambda,classes)
% Add L2 regularization to the objective function
[f, g] = logisticLoss(w,X,y,classes);
f = f + sum(lambda * (w.^2));
g = g + 2 * lambda * w;
end

function [yhat] = predict(model, Xhat)
% Predict labels using trained model
[nTest] = size(Xhat,1);

% Standardize optionally
if model.standardize
    Xhat = Xhat - repmat(model.mu, nTest, 1);
    Xhat = Xhat ./ repmat(model.sigma, nTest, 1);
end

if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end

% Predict labels (1..nClasses) for Xhat
num = Xhat * model.w;
[~, indx] = max(num, [], 2);
yhat = model.classes(indx);
end

function [lse] = logsumexp(b)
% Compute log(sum(exp)) across columns without overflowing
B = max(b,[],2);
lse = log(sum(exp(b-repmat(B,[1 size(b,2)])),2))+B;
end
