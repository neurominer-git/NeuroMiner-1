function [model] = ml_multiclass_bagging(X,y,options)
% ml_multiclass_bagging(X,y,options)
%
% Description:
%    - Classication based on the highest prediction among models trained on
%       bootstrap samples of dataset.
%
% Options:
%    - subModel: name of the underlying model (default: brokenStump)
%    - subOptions: options for subModel (default: none)
%    - nModels: number of models to fit the data to (default: 50)
%    - nSample: size of each boostrap sample (default: 0.7 * # of training data)
%
% Author:
%	 - Sampoorna Biswas (2014)

[nTrain, nFeatures] = size(X);

% Load default options
[subOptions, subModel, nModels, nSample] = ...
    myProcessOptions(options, 'subOptions', [], 'subModel', ...
    @ml_multiclass_stump, 'nModels', 50,'nSample', ...
    floor(nTrain * 0.7));

% Create models with random samples of trainig data
trainModels = cell(nModels, 1);
for k = 1:nModels
    train = randi([1 nTrain],nSample,1);
    XSampled{k} = X(train(:),:);
    ySampled{k} = y(train(:));
    trainModels{k} = subModel(X(train,:), y(train), subOptions);
end

% Set model fields
model.nModels = nModels;
model.predict = @predict;
model.XSampled = XSampled;
model.ySampled = ySampled;
model.trainModels = trainModels;
model.name = 'Classification with Bagging';
end

function [yhat] = predict(model, Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Predict and sum up individual models
yhats = zeros(nTest, model.nModels);
for k = 1:model.nModels
    training = model.trainModels{k};
    yhats(:,k) = training.predict(training, Xhat);
end

% Final result
yhat = mode(yhats, 2);
end
