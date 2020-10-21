function [model] = ml_binaryclass_stump(X,y,options)
% ml_binaryclass_stump(X,y,options)
%
% Description:
%    - Finds the best threshold across all features
%
% Options:
%	 - weights: define weights for every sample in training set (default: vector of 1's)
%
% Authors:
%	 - Nathaniel Lim (2014)

[nTrain, nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[z] = myProcessOptions(options, 'weights', ones(nTrain, 1));

% Select best feature and splitting threshold to minimize error
minErr = inf;
for j = 1:nFeatures
    sorted = sort(unique(X(:,j)));
    for t = 2:length(sorted)
        % Decide which side of the threshold belongs to which class by
        % minimizing error
        err1 = sum(z(y==-1).*(X(y==-1,j) < sorted(t))) + sum(z(y==1).*(X(y==1,j) >= sorted(t)));
        err2 = sum(z(y==1).*(X(y==1,j) < sorted(t))) + sum(z(y==-1).*(X(y==-1,j) >= sorted(t)));
        
        % Minimize error
        if err1 < minErr
            minErr = err1;
            minj = j;
            minThreshold = (sorted(t)+sorted(t-1))/2;
            thresholdType = 0;
        end
        if err2 < minErr
            minErr = err2;
            minj = j;
            minThreshold = (sorted(t)+sorted(t-1))/2;
            thresholdType = 1;
        end
    end
end

% Model outputs
model.name = 'Binary Decision Stump';
model.j = minj;
model.threshold = minThreshold;
model.thresholdType = thresholdType;
model.predict = @predict;
end

function yhat = predict(model, Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
yhat = zeros(nTest,1);

% Split data using threshold
index = Xhat(:,model.j) < model.threshold;

% Determine classes
if ~model.thresholdType
    yhat = sign(index - 0.5);
elseif model.thresholdType
    yhat = -sign(index - 0.5);
end

end