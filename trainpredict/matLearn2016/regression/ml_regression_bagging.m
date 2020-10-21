function model = ml_regression_bagging(X,y,options)
% ml_regression_bagging(X,y,options)
%
% Description:
%	 - Ensemble method that improves the stability and accuracy of regression algorithms. Also reduces variance 
%	   and helps to avoid overfitting.
%
% Options:
%    - subModel: regression model used for training (default: mean)
%    - subOptions: options used by subModel (default: none)
%    - nModels: number of bootstraps (default: 50)
%    - nSample: size of each bootstrap sample (default: 0.7 * # of training data)
%
% Authors:
%    - Jo�o Cardoso (2014)

[nTrain, nFeatures] = size(X);

% Load default options
[subOptions, subModel, nModels, nSample] = myProcessOptions(options, ...
    'subOptions', [], 'subModel', @ml_regression_mean, ...
    'nModels', 50,'nSample', floor(nTrain * 0.7));

% Create models with random samples of trainig data
trainModels = cell(nModels, 1);
for k = 1:nModels
    train = randi([1 nTrain],nSample,1);
    XSampled{k} = X(train,:);
    ySampled{k} = y(train);
    trainModels{k} = subModel(X(train,:), y(train), subOptions);
end

% Set model fields
model.name = 'Bagged Regression';
model.trainModels = trainModels;
model.nModels = nModels;
model.predict = @predict;
model.XSampled = XSampled;
model.ySampled = ySampled;
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

% Average results
yhat = ySum/model.nModels;
end
