function [model] = ml_binaryclass_perceptron(X,y,options)
% ml_binaryclass_perceptron(X,y,options)
%
% Description:
%	 - Fits a binary linear classifier by using the perceptron learning algorithm
%
% Options:
%	 - addBias: adds a bias variable (default: 1)
%    - lambdaL2: strength of L2-regularization parameter (default: 0)
%	 - alpha: subgradient descent learning rate (default: 1)
%	 - usePocket: uses the pocket algorithm to solve the stability problem
%       (default 0)
%	 - weights: gives the weights for each training example (default: vector of 1's)
%    - maxIter: maximum number of iterations; each iteration corresponds to
%       randomly selecting one misclassified example (default: 1e4)
%
% Authors:
% 	 - Kamyar Ardekani (2014)

[nTrain,nFeatures] = size(X);

% Set default options if not specified by input
[addBias,lambdaL2, alpha, usePocket, z, maxIter] = myProcessOptions(options,'addBias',1,'lambdaL2',0, 'alpha', 1, ...
    'usePocket', 0, 'weights', ones(nTrain, 1), 'maxIter', 1e4);

% Add optional bias
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

w0 = randn(nFeatures, 1);
w = w0;
minNumMisclassed = nTrain;
for iter = 1:maxIter
    % Prediction
    index = X*w >= 0;
    yhat = sign(index - 0.5);
    
    misclassedIndices = yhat ~= y;
    % Find misclassified points
    XMisclassed = X(misclassedIndices, :);
    yhatMisclassed = yhat(misclassedIndices);
    nMisclassed = length(yhatMisclassed);
    
    % Look for the w that gives the minimum number of misclassified
    % examples (for pocket learning algorithm)
    if usePocket
        if nMisclassed < minNumMisclassed
            minNumMisclassed = nMisclassed;
            bestW = w;
        end
    end
    
    % Stop if all data points are classified correctly
    if isempty(XMisclassed)
        break
    end
    
    % Randomly select a misclassified example to update w based on that example
    updateIndex = randsample(nMisclassed, 1, true, z(misclassedIndices));
    
    % Update w using subgradient descent method
    % Bias parameter w(0) is not penalized
    w = w - alpha * (yhatMisclassed(updateIndex) * XMisclassed(updateIndex, :)' + [0; lambdaL2*w(2:end)]);
end

% Choose best weights with pocket algorithm
if usePocket
    w = bestW;
end

% Model outputs
model.name = 'Perceptron Binary Classification';
model.addBias = addBias;
model.w = w;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
yhat = size(nTest,1);

% Add optional bias
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end

% Final results
index = Xhat*model.w >= 0;
yhat = sign(index - 0.5);
end
