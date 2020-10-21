function [model] = ml_multiclass_boosting(X,y,options)
% ml_multiclass_boosting(X,y,options)
%
% Description:
%	 - Implements a SAMME AdaBoosting algorithm to boost a classifier 
%      (see in Zhu, Rosset, Zou, and Hastie 2006)
%
% Options:
%    - nBoosts: the number of weak learners to train (default: 10)
%    - subModel: the classifier to use (default: decision stump)
%    - subOptions: options (other than weights) passed into subModel
%
% Authors:
% 	 - Rebecca McKnight (2014)

[nTrain,nFeatures] = size(X);
classes = unique(y);
nClasses = length(classes);

% Process options
[nBoosts,subModel, subOptions] = myProcessOptions(options, 'nBoosts',50, ...
    'subModel',@ml_multiclass_stump, 'subOptions', []);

% Initialize weights
z = (1/nTrain)*ones(nTrain,1);
alpha = zeros(nBoosts,1);

for k = 1:nBoosts
    % Train weighted classifier
    subOptions.weights = z;
    model.trainModels{k} = subModel(X,y,subOptions);
    
    % Compute predictions
    yhat = model.trainModels{k}.predict(model.trainModels{k},X);
    
    % Compute weighted error rate
    err = (z'*(y~=yhat))/sum(z);
    
    % Compute alpha
    alpha(k) = log((1-err)/err) + log(nClasses - 1);
    
    % Update weights
    z = z.*exp(alpha(k).*(y~=yhat));
    
    % Re-normalize weights
    z = z/sum(z);
end


model.name = 'AdaBoosted Classification';
model.alpha = alpha;
model.nBoosts = nBoosts;
model.nClasses = nClasses;
model.classes = classes;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);
alpha = model.alpha;

% Predict and sum up individual models
yhats = zeros(nTest, model.nClasses);
for k = 1:model.nBoosts
    training = model.trainModels{k};
    pred = training.predict(training, Xhat);
    yBin = bsxfun(@eq, pred, model.classes');
    yhats = yhats + alpha(k)*yBin;
end

% Final result
[~, indx] = max(yhats, [], 2);
yhat = model.classes(indx);
end
