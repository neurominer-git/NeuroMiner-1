function [ model ] = ml_regression_local(X,y,options)
% matlearn_regression_local(X,y,options)
%
% Description:
%    - Performs local regression, a meta-model based on k-nearest-neighbor,
%       where neighbouss are weighted according to distances & using a subModel
%       trained on subset
%    - Requires a fairly large and dense data to produce a good
%       model; does not produce a function easily represented by a
%       mathematical formula; and is quite computationally expensive
%
% Options:
%	 - k: number of points to consider in each local data set
%    - subModel: regression model to apply in each local data set (default: mean)
%    - subOptions: options to pass to the subModel; note that
%       local regression will overwrite the 'weights' option (default: none)
%    - weightingFunc: function usd to calculate weights for the local data set,
%       takes in vector of distances between 0 and 1, of length k
%       (default: @triCube)
%
% Authors:
%    - Giorgio Gori (2014)


% Get and save options
[subModel, subOptions, k, weightingFunc] = ...
    myProcessOptions(options, 'subModel', @ml_regression_mean, ...
    'subOptions', [], 'k', 1, 'weightingFunc', @ml_weights_triCube);


% Save the training set
model.X = X;
model.y = y;

% Options
model.subModel = subModel;
model.subOptions = subOptions;
model.k = k;
model.weightingFunc = weightingFunc;

% Name and function definition
model.name = 'Local Regression';
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
% Recover the training data and options from the model
X = model.X;
y = model.y;
k = model.k;

[nTrain, nFeatures] = size(X);
[nTest,nFeatures] = size(Xhat);

% Pre-compute all the distances
D = X.^2*ones(nFeatures,nTest) + ones(nTrain,nFeatures)*(Xhat').^2 - 2*X*Xhat';

yhat = zeros(nTest,1);
for i = 1:nTest
    % Find the closest k neighbors
    [minDist,sorted] = sort(D(:,i));
    neighbors = sorted(1:k);
    
    % If the test point is the same as a train point don't divide by 0
    if minDist(k) == 0
        minDist(k) = 1;
    end
    
    % Compute the weights using weightingFunc
    model.subOptions.weights = model.weightingFunc(minDist(1:k)/minDist(k));
    
    % If the sum of weights is zero, make them ones
    if sum(model.subOptions.weights) == 0
        model.subOptions.weights = ones(k, 1);
    end
    
    % Fit sub model on this local data set
    local_subModel = model.subModel(X(neighbors, :), y(neighbors, :), model.subOptions);
    yhat(i) = local_subModel.predict(local_subModel, Xhat(i,:));
end
end
