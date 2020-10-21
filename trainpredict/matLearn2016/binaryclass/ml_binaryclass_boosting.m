function [model] = ml_binaryclass_boosting(X,y,options)
% ml_binaryclass_boosting(X,y,options)
%
% Description:
%	 - Implements a boosting algorithm to boost a binary classifier
%
% Options:
%    - nBoosts: the number of weak learners to train (default: 10)
%    - subModel: the classifier to use (default: decision stump)
%    - booster: the boosting algorithm to use, one of:
%       'ada' - AdaBoost
%       'logit' - LogitBoost
%       (default: ada)
%    - subOptions: options (other than weights) passed into subModel
%
% Authors:
% 	- Rebecca McKnight (2014)

[nTrain,nFeatures] = size(X);

% Process options
[nBoosts,subModel,booster, subOptions] = myProcessOptions(options, 'nBoosts',50, ...
    'subModel',@ml_binaryclass_stump, ...
    'booster','ada', 'subOptions', []);

% Initialize weights
z = (1/nTrain)*ones(nTrain,1);
alpha = zeros(nBoosts,1);

for k = 1:nBoosts
    % Train weighted classifier
    subOptions.weights = z;
    model.trainModels{k} = subModel(X,y,subOptions);
    
    % Compute predictions
    p = model.trainModels{k}.predict(model.trainModels{k},X);
    
    % Compute weighted error rate
    if isstruct(p), yhat= p.yhat; else yhat = p; end
    err = sum(z.*(y~=yhat));
    
    % Compute alpha
    alpha(k) = (1/2)*log((1-err)/err);
    
    % Update weights
    if strcmp(booster,'logit')
        % Logitboost
        z = z.*log(ones(nTrain,1) + exp(-alpha(k)*y.*yhat));
    else
        % Default to Adaboost
        z = z.*exp(-alpha(k)*y.*yhat);
    end
    
    % Re-normalize weights
    z = z/sum(z);
end

if strcmp(booster,'logit')
    model.name = 'LogitBoosted Binary Classification';
else
    model.name = 'AdaBoosted Binary Classification';
end
model.alpha = alpha;
model.nBoosts = nBoosts;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);
alpha = model.alpha;

% Predict and sum up individual models
ySum = zeros(nTest, 1);
for k = 1:model.nBoosts
    training = model.trainModels{k};
    p = training.predict(training,Xhat);
    if isstruct(p), yhat= p.yhat; else yhat = p; end
    ySum = ySum + alpha(k)*yhat;
end

% Final result
index = ySum >= 0;
yhat = sign(index - 0.5);
end