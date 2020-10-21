function [model] = ml_multiclass_decisionTree(X,y,options)
% ml_multiclass_decisionTree(X,y,options)
%
% Description:
%       - A decision tree that classifies data into multiple class labels
% Options:
%   - maxDepth: The max depth of the decision tree (default: 8)
%   - minLeafSize: The min lea node size of the decision tree (default: 3)
%   - splitType: how to choose the feature and threshold to split on
%       'info' - information gain
%       'err' - classification error
%       (default: 'info')
%
% Authors:
% 	- Jeff Allen (2014)

% Load options
[maxDepth,minLeafSize, splitType] = myProcessOptions(options,...
    'maxDepth', 8, ...
    'minLeafSize', 3, ...
    'splitType','info');

% Unique classes of y
classes = unique(y);
options.maxDepth = maxDepth;
options.minLeafSize = minLeafSize;
options.splitType = splitType;
options.classes = classes;

% Recursively split data and build the tree
root = initNode(X, y, 0);
root = buildTree(root,options);

% Model ouputs
model = [];
model.name = 'Decision Tree';
model.root = root;
model.predict = @predict;
end

function [node] = buildTree(node,options)
X = node.X;
y = node.y;
[nTrain, nFeatures] = size(X);

% Return most common label
if node.depth >= options.maxDepth || nTrain <= options.minLeafSize ||length(unique(y)) == 1;
    node.val = mode(y);
    node.isLeaf = 1;
    return;
else
    % Determine feature with largest information gain and create child
    % nodes
    [node] = optimalSplit(node, options);
    leftChild = buildTree(node.leftChild, options);
    rightChild = buildTree(node.rightChild, options);
    node.leftChild = leftChild;
    node.rightChild = rightChild;
end
end
function [node] = optimalSplit(node, options)
% Fits the tree based on the options provided
X = node.X;
y = node.y;

[nTrain,nFeatures] = size(X);

if strcmp(options.splitType,'info')
    % Split based on information gain
    maxInfoGain = -inf;
    for j = 1:nFeatures
        sorted = sort(unique(X(:,j)));
        for t = 2:length(sorted)
            infoGain = calcInfoGain(j,sorted(t),X,y,options.classes);
            if infoGain > maxInfoGain
                maxInfoGain = infoGain;
                optimalThreshold = (sorted(t) + sorted(t-1))/2;
                optimalDimension = j;
            end
        end
    end
elseif strcmp(options.splitType, 'err')
    % Split based on classification error
    minErr = inf;
    for j = 1:nFeatures
        sorted = sort(unique(X(:,j)));
        for t = 2:length(sorted)
            err1 = sum(X(y==mode(y),j) < sorted(t));
            err2 = sum(X(y~=mode(y),j) >= sorted(t));
            err = err1 + err2;
            
            if err < minErr
                minErr = err;
                optimalDimension = j;
                optimalThreshold = (sorted(t) + sorted(t-1))/2;
            end
        end
    end
end

index = (X(:,optimalDimension) < optimalThreshold);
leftX  = X(index, :);
lefty  = y(index , :);
rightX = X(~index, :);
righty = y(~index, :);

node.leftChild = initNode(leftX,lefty,node.depth + 1);
node.rightChild = initNode(rightX,righty,node.depth + 1);

node.splitDimension = optimalDimension;
node.threshold = optimalThreshold;
end

function [infoGain] = calcInfoGain(featureDim,threshold,X,y,classes)
% Returns the info gain of the data set given a feature to split
% and a threshold to test

[nTrain,nFeatures] = size(X);

% Entropy of original set
totalEntropy = calcEntropy(y,classes);

% Get all the indicies less than the threshold for that feature
index = X(:,featureDim) < threshold;

lefty = y(index);
nLeft = sum(index);
righty = y(~index);
nRight = sum(~index);

% Entropy of the two subsets
leftEntropy = calcEntropy(lefty,classes);
rightEntropy = calcEntropy(righty,classes);

infoGain = totalEntropy-(nLeft/nTrain)*leftEntropy-(nRight/nTrain)*rightEntropy;
end

function [entropy] = calcEntropy(y,classes)
% Calculate the entropy of a data set
[classProbs] = getClassProbs(y,classes);
logProbs = zeros(length(classes),1);
logProbs(classProbs~= 0) = log2(classProbs(classProbs~= 0));

entropy = -sum(classProbs'*logProbs);
end

function [classProbs] = getClassProbs(y,classes)
% Returns the class probabilities (# occurrences / total training data)
[nTrain,nFeatures] = size(y);
nClasses = length(classes);

% Count each class's occurrence
classCounts = zeros(nClasses,1);
for c = 1 : nClasses
    classCounts(c) = sum(y==classes(c));
end
classProbs = classCounts/nTrain;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Recurse until leaf nodes to detemine yhat
yhat = predictTree(Xhat,-inf*ones(nTest, 1),model.root, 1:nTest);
end

function [yhat] = predictTree(Xhat,yhat,node,indices)
% Predicts the classes based on the decision tree
[nTest,nFeatures] = size(Xhat);

if node.isLeaf == 1
    % Determine yhat if node is leaf node
    if length(indices)>0
        yhat(indices) = node.val;
    end
else
    childInd1 = intersect(indices,find(Xhat(:,node.splitDimension) < node.threshold));
    childInd2 = intersect(indices,find(Xhat(:,node.splitDimension) >= node.threshold));
    
    childyhat1 = predictTree(Xhat,yhat,node.leftChild,childInd1);
    childyhat2 = predictTree(Xhat,yhat,node.rightChild,childInd2);
    yhat = max(childyhat1, childyhat2);
end
end

function [node] = initNode(X, y, depth)
% Initialize node
node = [];
node.X = X;
node.y = y;

% Common tree structure
node.leftChild = [];
node.rightChild = [];
node.depth = depth;
node.isLeaf = 0;
node.splitDimension = 1;
node.val = 0;
node.type = 0;
node.threshold = inf;
end
