function [model] = ml_regression_KNN(X,y,options)
% ml_regression_KNN(X,y,options)
%
% Description:
%    - Interpolation based k-nearest neighbours (default: 10)
%
% Options:
%    - k: number of nearest neighbours
%
% Authors: Fujun Xie (2014)

% Use default value of options if not defined
[k] = myProcessOptions(options,'k', 10);

% Set model outputs
model.x = X;
model.y = y;
model.name = 'K-Nearest Neighbour Regression';
model.predict = @predict;
model.k = k;
end

function [yhat] = predict(model,Xhat)
X = model.x;
y = model.y;
[nTrain,nFeatures] = size(X);
[nTest,nFeatures] = size(Xhat);
% Calculate all Euclidean distances between X and Xhat
D = X.^2*ones(nFeatures,nTest) + ones(nTrain,nFeatures)*(Xhat').^2 - 2*X*Xhat';
yhat = zeros(nTest,1);
% Calculate mean of k nearest neighbours for each test data
for i = 1:nTest
    [dist,sorted] = sort(D(:,i));
    yhat(i) = sum(y(sorted(1:model.k)))/model.k;
end
end