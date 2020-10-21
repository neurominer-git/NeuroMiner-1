function [model] = ml_multiclass_KNN(X,y,options)
% ml_multiclass_KNN(X,y,options)
%
% Description:
%	 - An object is classified by a majority vote of its neighbors, with the object being assigned
%       to the class most common among its k nearest neighbors.
%
% Options:
%    - k : number of nearest neighbors to consider for finding
%   the class of each example (default: 10)
%
% Authors:
%    - Nasim (Sedigheh) Zolaktaf (2014)

% Sets the default values which are not set in the option
[k] = myProcessOptions(options,'k', 10);

% Model outputs
model.name = 'k-Nearest Neighbours Classification';
model.k = k;
model.X = X;
model.y = y;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% returns yhat, where yhat(i) is equal to the class training example i belongs to
X= model.X;
y = model.y;
[nTrain, nFeatures] = size(X);
[nTest, nFeatures] = size(Xhat);

% Calculate all Euclidean distances between X and Xhat
D = X.^2*ones(nFeatures,nTest) + ones(nTrain,nFeatures)*(Xhat').^2 - 2*X*Xhat';

yhat = zeros(nTest,1);
% Calculate mode of k nearest neighbours for each test data
for i = 1:nTest
    [sortedValues,sortIndex] = sort(D(:,i));
    yhat(i) = mode(y(sortIndex(1:model.k)));
end
end
