function [model] = ml_binaryclass_bagging(X,y,options)
% ml_binaryclass_bagging(X,y,options)
%
% Description:
%    - Classification based on the highest prediction among models to bootstrap samples.
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
[subOptions, subModel, nModels, nSample] = myProcessOptions(options,...
    'subOptions', [], 'subModel', @ml_binaryclass_brokenStump,...
    'nModels', 50,'nSample', floor(nTrain * 0.7));

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
model.name = sprintf('Binary Classification with Bagged %s',...
                     trainModels{1}.name);
end

function [yhat] = predict(model, Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Predict and sum up individual models
ySum = zeros(nTest, 1);
for k = 1:model.nModels
    training = model.trainModels{k};
    ySum = ySum + training.predict(training, Xhat);
end

% Final result
index = ySum >= 0;
yhat = sign(index - 0.5);
end
