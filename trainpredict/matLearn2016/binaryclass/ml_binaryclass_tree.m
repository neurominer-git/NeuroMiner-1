function [model] = ml_binaryclass_tree(X,y,options)
% ml_binaryclass_tree(X,y,options)
%
% Describtion:
%	 - Binary classification using a decision tree. The tree is constructed
%       using the ID3/C4.5 algorithm, based on information gain
%
% Options:
%    - dataType: vector 0/1's where 0 indicate i'th feature is discrete,
%       and 1 means feature is continuous (default: vector of 1's)
%    - maxDepth: The maximum depth of the tree (default: 8)
%    - minLeafSize: The mininum number of data in the leaf node
%       (default: 3)
%
% Authors:
%    - Philipp Witte (2014)

[nTrain, nFeatures] = size(X);

% Check optional parameters and set default parameter
[maxDepth,minLeafSize,dataTypes] = myProcessOptions(options,'maxDepth', 8, 'minLeafSize', 3,'dataTypes', ones(nFeatures,1));

options.maxDepth = maxDepth;
options.minLeafSize = minLeafSize;
options.dataTypes = dataTypes;

% Recursively split data and build the tree
root = initNode(X,y,0);
root = buildTree(root,options,zeros(nFeatures,1));

% Output model
model = [];
model.name = 'Binary Decision Tree';
model.root = root;
model.dataTypes = options.dataTypes;
model.predict = @predict;
end

function [node] = buildTree(node,options,isUsed)
X = node.X;
y = node.y;
[nTrain, nFeatures] = size(X);

% Return most common label
if all(isUsed) || node.depth >= options.maxDepth || nTrain <= options.minLeafSize ||length(unique(y)) == 1;
    node.val = mode(y);
    node.isLeaf = 1;
    return;
else
    % Determine feature with largest information gain and create child
    % nodes
    [node, isUsed] = optimalSplit(node, options, isUsed);
    for i = 1:length(node.children)
        child = buildTree(node.children{i}, options, isUsed);
        node.children{i} = child;
    end
end
end

function [node, isUsed] = optimalSplit(node,options,isUsed)
X = node.X;
y = node.y;

% Determine feature for new tree branch using information gain
[nTrain,nFeatures] = size(X);

maxInfoGain = -inf;
for j=1:nFeatures
    if ~isUsed(j) && options.dataTypes(j) == 0;
        infoGain = calcInfoGainDisc(X(:,j),y);
        if infoGain > maxInfoGain;
            maxInfoGain = infoGain;
            optimalDimension = j;
        end
    elseif options.dataTypes(j) == 1
        sorted = sort(unique(X(:,j)));
        for t = 2:length(sorted)
            infoGain = calcInfoGainCont(sorted(t),X(:,j),y);
            if infoGain > maxInfoGain;
                maxInfoGain = infoGain;
                optimalThreshold = (sorted(t)+sorted(t-1))/2;
                optimalDimension = j;
            end
        end
    end
end

% Create child nodes based on optimal feature
if (options.dataTypes(optimalDimension) == 0)
    uniqueVals = unique(X(:,optimalDimenions));
    for v=1:length(uniqueVals)
        index = X(:,optimalDimension) == uniqueVals(v);
        child = initNode(X(index,:),y(index,:),node.depth + 1);
        child.class = uniqueVals(v);
        node.children{v} = child;
    end
    isUsed(optimalDimension) = 1;
    node.type = 0;
elseif(options.dataTypes(optimalDimension) == 1)
    index = (X(:,optimalDimension) < optimalThreshold);
    leftX  = X(index, :);
    lefty  = y(index , :);
    rightX = X(~index, :);
    righty = y(~index, :);
    
    leftChild = initNode(leftX,lefty,node.depth + 1);
    rightChild = initNode(rightX,righty,node.depth + 1);
    node.children{1} = leftChild;
    node.children{2} = rightChild;
    node.threshold = optimalThreshold;
    node.type = 1;
end
node.splitDimension = optimalDimension;
end

function [infoGain] = calcInfoGainCont(threshold,X,y)
% Returns the info gain of the data set given a feature to split
% and a threshold to test
[nTrain,~] = size(X);
classes = unique(y);
countLabel1 = sum(y(:) == classes(1));
probClass1 = countLabel1/nTrain;
probClass2 = 1 - probClass1;

% Entropy of original set
totalEntropy = -probClass1*log2(probClass1) - probClass2*log2(probClass2);

index = X < threshold;
nLeft = sum(index);
nLeftClass1 = sum(index & (y(:) == classes(1)));
nLeftClass2 = nLeft - nLeftClass1;
nRight = sum(~index);
nRightClass1 = sum(~index & (y(:) == classes(1)));
nRightClass2 = nRight - nRightClass1;

if ((nLeftClass1 == 0) || (nLeftClass2 == 0))
    % Feature of value v already classified
    leftEntropy = 0;
else
    % Entropy of subset
    leftEntropy = -nLeftClass1/nLeft*log2(nLeftClass1/nLeft) - ...
        nLeftClass2/nLeft*log2(nLeftClass2/nLeft);
end
if ((nRightClass1 == 0) || (nRightClass2 == 0))
    % Feature of value v already classified
    rightEntropy = 0;
else
    % Entropy of subset
    rightEntropy = -nRightClass1/nRight*log2(nRightClass1/nRight) - ...
        nRightClass2/nRight*log2(nRightClass2/nRight);
end
infoGain = totalEntropy - (nLeft/nTrain)*leftEntropy - ...
                 (nRight/nTrain)*rightEntropy;
end

function [infoGain] = calcInfoGainDisc(X,y)
% Returns the info gain of the data set given a feature to split
% and a threshold to test
[nTrain,~] = size(X);
classes = unique(y);
countLabel1 = sum(y(:) == classes(1));
probClass1 = countLabel1/nTrain;
probClass2 = 1 - probClass1;

totalEntropy = -probClass1*log2(probClass1) - probClass2*log2(probClass2);
infoGain = totalEntropy;

uniqueVals = unique(X);
for v=1:length(uniqueVals)
    index = (X == uniqueVals(v));
    nSubset = sum(index);
    nClass1 = sum(index & (y(:) == classes(1)));
    nClass2 = nSubset - nClass1;
    if ((nClass1 == 0) || (nClass2 == 0))
        % Feature of value v already classified
        subsetEntropy = 0;
    else
        % Entropy of subset
        subsetEntropy = -nClass1/nSubset*log2(nClass1/nSubset) - ...
            nClass2/nSubset*log2(nClass2/nSubset);
    end
    infoGain = infoGain - (nSubset/nTrain)*subsetEntropy;
end
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Recurse until leaf nodes to detemine yhat
yhat = predictTree(Xhat,-inf*ones(nTest, 1),model.root, 1:nTest);
end

function [yhat] = predictTree(Xhat,yhat,node,indices)
% Recursively travel from root to leaf to predict values

if node.isLeaf == 1
    % Determine yhat if node is leaf node
    yhat(indices) = node.val;
else
    % Set yhat to values changed by child nodes
    if node.type == 0
        for i = 1:length(node.children)
            childInd = union(indices, find(Xhat(:,node.splitDimension) == node.children{i}.class));
            childyhat = predictTree(Xhat,yhat,node.children{i},childInd);
            yhat = max(yhat, childyhat);
        end
    elseif node.type == 1
        childInd1 = intersect(indices,find(Xhat(:,node.splitDimension) < node.threshold));
        childInd2 = intersect(indices,find(Xhat(:,node.splitDimension) >= node.threshold));
        
        childyhat1 = predictTree(Xhat,yhat,node.children{1},childInd1);
        childyhat2 = predictTree(Xhat,yhat,node.children{2},childInd2);
        yhat = max(childyhat1, childyhat2);
    end
end
end

function [node] = initNode(X, y, depth)
% Initialize node
node = [];
node.X = X;
node.y = y;

% Common tree structure
node.children = [];
node.depth = depth;
node.isLeaf = 0;
node.splitDimension = 1;
node.val = 0;
node.type = 0;

% For real features
node.threshold = inf;

% For finite discrete features
node.class = 0;
end