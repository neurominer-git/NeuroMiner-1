function [model] = ml_binaryclass_randomForest(X,y,options)
% ml_binaryclass_randomForest(X,y,options)
%
% Description:
%    - A bagging of tree classifers, where the splits are chosen
%        among random selection of features
% Options:
%    - nModels: number of trees in the random forest (default: 10)
%    - nSample: size of each boostrap sample (default: 0.7 * # of training data)
%    - minLeafSize: the minimum number of training data in leaf nodes
%       (default: 3)
%    - maxDepth: the max depth of each decision tree (default: 8)
%    - maxFeatures: the number of features that each decision tree will
%                 consider (default: sqrt(nFeatures))
%    - oobMaxFeatures: a list of maxFeatures to try; the best one
%       chosen by out of bag error
%       * overides the maxFeatures input
%
% Output:
%    - oobError: out of bag error. Estimates genealized error in a similar
%       way to cross validation does
%    - oobPred: final out of bag predictions
%    - oob: matrix of all OOB predictions for each tree (0 if used for
%       training)
%
% Authors:
%  	 - Alim Virani (2014)

[nTrain, nFeatures] = size(X);

if isfield(options,'nSample'), options.nSample = floor(nTrain * options.nSample); end
if isfield(options,'maxFeatures'), options.maxFeatures = floor(nFeatures * options.maxFeatures); end
% Load default options
[nModels, maxFeatures, minLeafSize, maxDepth, oobMaxFeatures, nSample] = myProcessOptions(options,'nModels',10,'maxFeatures',sqrt(nFeatures), 'minLeafSize', 3, 'maxDepth', 8, 'oobMaxFeatures', [], 'nSample', floor(nTrain * 0.7));
treeOptions.maxDepth = maxDepth;
treeOptions.minLeafSize = minLeafSize;

if isempty(oobMaxFeatures)
    % Build a forest with maxFeatures
    oob = zeros(nTrain, nModels);
    for k = 1:nModels
        
        % Create models with random samples of trainig data
        train = randi([1 nTrain],1,nSample);
        XSampled{k} = X(train,:);
        ySampled{k} = y(train);
        featInd = randperm(nFeatures);
        features{k} = featInd(1:maxFeatures);
        
        trainModel = ml_binaryclass_tree(X(train,features{k}),y(train),treeOptions);
        
        % OOB predictions for tree k
        % OOB predictions are done inefficiently and implementation should be
        % improved
        yhat = trainModel.predict(trainModel, X(:,features{k}));
        yhat(train) = 0;
        oob(:,k) = yhat;
        
        trainModels{k} = trainModel;
    end
    
    % Out of bag error for cross-validation Of maxFeatures
    index = sum(oob,2) >= 0;
    oobPred = sign(index - 0.5);
    oobError = sum(y~=oobPred);
    
    model.predict = @predict;
    model.oob = oob;
    model.oobPred = oobPred;
    model.oobError = oobError;
    model.name = 'Random Forest Binary Classification';
    model.features = features;
    model.XSampled = XSampled;
    model.ySampled = ySampled;
    model.trainModels = trainModels;
    model.nModels = nModels;
else
    % Build a RandomForest for each maxFeatures value in OOBmaxFeatures
    minOOBError = inf;
    for i = 1:length(oobMaxFeatures)
        optionsCandidate.nModels = nModels;
        optionsCandidate.maxDepth = maxDepth;
        optionsCandidate.minLeafSize = minLeafSize;
        optionsCandidate.maxFeatures =oobMaxFeatures(i);
        
        candidateModel = ml_binaryclass_randomForest(X,y,optionsCandidate);
        
        if candidateModel.oobError<minOOBError
            % Current best error, and associated model
            minOOBError = candidateModel.oobError;
            model = candidateModel;
        end
    end
end
end

function p = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);

% Predict and sum up individual models
p.D = zeros(nTest, 1);
for k = 1:model.nModels
    training = model.trainModels{k};
    p.D = p.D + training.predict(training, Xhat(:,model.features{k}));
end

% Final result
index = p.D >= 0;
p.yhat = sign(index - 0.5);
end
