function [model] = ml_regression_tree(X,y,options)
% ml_regression_tree(X,y,options)
%
% Description:
%    - Minimizes the squared error using a regression tree
%
% Options:
%    - maxDepth: The maximum depth of the tree (default: 8)
%    - minLeafSize: The mininum number of data in the leaf node
%       (default: 8)
%	 - modelType: 'cns' or 'lin' model in leaf node
%       (default: 'lin')
%    - lambdaL2: L2 regularization weigting for 'linear' modelType
%       (default: 0)
%    - splitSample: method of splitting the samples, one of:
%       'bf' - brute force, minimizing squared error of difference between
%       subset values and mean
%       'rnd' - randomly selects thresholds for splitting data, much faster
%       (default: 'rnd')
%
% Authors:
% 	 - Jianhui Chen (2014)

% Check optional parameters and set default parameter
[maxDepth, minLeafSize, lambdaL2, modelType, splitSample] = myProcessOptions(options, 'maxDepth', 8, 'minLeafSize', 8, 'lambdaL2', 0, 'modelType', 'lin', 'splitSample', 'rnd');

options.maxDepth = maxDepth;
options.minLeafSize = minLeafSize;
options.lambdaL2 = lambdaL2;
options.splitSample = splitSample;
options.modelType = modelType;

% Recursively split data and build the tree
root = initNode(X, y, 0);
root = buildTree(root, options);

% Output model
model = [];
model.name = 'Regression Tree';
model.predict = @predict;
model.root = root;
model.modelType = modelType;
end

function [node] = buildTree(node, options)
% Recursive function to build a tree

% Stop splitting when tree is too deep or the data size too small
if node.depth >= options.maxDepth || length(node.y) <= options.minLeafSize
    node.isLeaf = 1;
    if strcmp(options.modelType, 'cns')
        node.mean = mean(node.y);
    elseif strcmp(options.modelType, 'lin')
        XwithBias = [ones(size(node.X, 1), 1), node.X];
        [w, yhat] = leastSquareL2normWeight(XwithBias, node.y, options.lambdaL2);
        node.w = w;
    end
    return;
end

% Otherwise split current node
[node] = optimalSplit(node, options);

[leftChild] = buildTree(node.leftChild, options);
node.leftChild = leftChild;

[rightChild] = buildTree(node.rightChild, options);
node.rightChild = rightChild;

end

function [node] = optimalSplit(node, options)
% Split and initialize node's leftChild and rightChild
X = node.X;
y = node.y;
[nTrain, nFeatures] = size(X);

minError = inf;
optimalDimension = 1;
optimalThreshold = 0;
% Loop over all dimensions to pick the optimal dimension to split data
for i=1:nFeatures
    [curThreshold, curError] = minimizeRegressionError(X(:,i), y, nFeatures, options);
    
    if curError < minError
        minError = curError;
        optimalDimension = i;
        optimalThreshold = curThreshold;
    end
end

% Split data by given dimension and threshold
index = (X(:, optimalDimension) <= optimalThreshold);
leftX  = X(index, :);
lefty  = y(index , :);
rightX = X(~index, :);
righty = y(~index, :);

% Create leftChild node and rightChild node
leftChild  = initNode(leftX, lefty, node.depth + 1);
rightChild = initNode(rightX, righty, node.depth + 1);
node.splitDimension = optimalDimension;
node.threshold = optimalThreshold;
node.leftChild  = leftChild;
node.rightChild = rightChild;
end

function [threshold, minErr] = minimizeRegressionError(X, y, featureDim, options)
% Optimal threshold to minimize (squared) regression error
[nTrain, ~] = size(X);
if strcmp(options.splitSample, 'bf')
    % Brute force sample possible splitting values
    % Complexity O(Nlog(N)),  N is data size for whole tree
    sorted = sort(unique(X));
    % Loop over all possible splitting values
    minErr = inf;
    for t = 2:length(sorted)
        index = (X < sorted(t));
        yLeft  = y(index);
        yRight = y(~index);
        
        % Squared fitting error in left and right child, choose minimal error
        err = sum((yLeft - mean(yLeft)).^2) + sum((yRight - mean(yRight)).^2);
        if err < minErr
            threshold = (sorted(t)+sorted(t-1))/2;
            minErr = err;
        end
    end
elseif strcmp(options.splitSample, 'rnd')
    % Randomly sample # featureDim*2 value in the range [x_min x_max]
    % Complexity O(dlog(N)) d is feature dimension, N is data size for whole tree
    XMin = min(X);
    XMax = max(X);
    thresholds = rand(1, featureDim*2)*(XMax - XMin) + XMin;
    % Loop over all possible splitting values
    minErr = inf;
    for t = thresholds
        index = (X < t);
        yLeft  = y(index);
        yRight = y(~index);
        % Squared fitting error in left and right child, choose minimal error
        err = sum((yLeft - mean(yLeft)).^2) + sum((yRight - mean(yRight)).^2);
        if err < minErr
            threshold = t;
            minErr = err;
        end
    end
end
end

function [w, yhat] = leastSquareL2normWeight(X, y, lambda)
% Creates linear regression (L2) model with training sample weights and L2
% regularization

% Optimization options for minFunc
optimOptions.Display = 0;
optimOptions.useMex = 0;

% compute weight vector w using minFunc
w = minFunc(@squaredLossL2,randn(size(X,2),1),optimOptions,X,y,lambda);
yhat = X * w;
end

function [f,g] = squaredLossL2(w,X,y,lambda)
% Squared loss function with training sample weights and L2 regularization
f = (X*w-y)'*(X*w-y) + lambda*(w'*w);
g = 2*X'*(X*w-y) + 2*lambda*w;
end

function [yhat] = predict(model, Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);

% Recurse until leaf nodes to detemine yhat
yhat = predictTree(Xhat, -inf*ones(nTest,1), model.root, 1:nTest, model.modelType);
end

function [yhat] = predictTree( Xhat, yhat, node, indices, modelType)
% Recursively travel from root to leaf to predict values

if node.isLeaf
    % Determine yhat if node is leaf node
    if strcmp(modelType, 'cns')
        % Constant model
        yhat(indices) = node.mean;
    elseif strcmp(modelType, 'lin')
        % Linear model
        Xhat = [ones(size(Xhat,1),1), Xhat];
        yhat(indices) = Xhat(indices,:)*node.w;
    end
    return;
else
    % Set yhat to values changed by child nodes
    childIndLeft = intersect(indices,find(Xhat(:,node.splitDimension) < node.threshold));
    childIndRight = setdiff(indices, childIndLeft);
    childyhatLeft = predictTree(Xhat, yhat, node.leftChild, childIndLeft, modelType);
    childyhatRight = predictTree(Xhat, yhat, node.rightChild, childIndRight, modelType);
    yhat = max(childyhatLeft, childyhatRight);
    
end
end

function [node] = initNode(X, y, depth)
% Initialize node
node.X = X;
node.y = y;

% Common tree structure
node.leftChild = struct();
node.rightChild = struct();
node.depth = depth;
node.isLeaf = 0;
node.splitDimension = 1;
node.threshold = 0;

% Constant leaf node
node.mean = 0;
% Linear leaf node
node.w = 0;
end
