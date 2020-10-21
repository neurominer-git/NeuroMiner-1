function [model] = ml_multiclass_stump(X,y,options)
% ml_multiclass_stump(X,y,options)
%
% Description:
%	 - Finds the best threshold across all features
%
% Options:
%    - weights: weights for training data (default: vector of 1's)
%
% Authors:
% 	 - Antoine Ponsard (2014)

[nTrain,nFeatures] = size(X);

% Input options
[z] = myProcessOptions(options, 'weights', ones(nTrain, 1));

% Number of unique classes
classes = unique(y);
nClasses = length(classes);

% number of training data points that fall in each leaf
% For a stump, only two leaves: one below and one above the threshold
count = zeros(nClasses, 2);

% Error to be minimized at each step
minErr = inf;
for j = 1:nFeatures
    % Try every possible threshold based on the values of X along feature j
    sorted = sort(unique(X(:,j)));
    for t = 2:length(sorted)
        for k=1:nClasses
            count(k, 1) = sum(y(X(:,j)<sorted(t)) == classes(k));
            count(k, 2) = sum(y(X(:,j)>=sorted(t)) == classes(k));
        end
        
        count(:,1) = count(:,1)./(sum(count(:,1)));
        count(:,2) = count(:,2)./(sum(count(:,2)));
        
        % Most likely class for each leaf and its empirical probability
        [pyhat(1), ind1] = max(count(:,1));
        [pyhat(2), ind2] = max(count(:,2));
        yhat(1) = classes(ind1);
        yhat(2) = classes(ind2);
        err1 = sum(z(y~=yhat(1)).*(X(y~=yhat(1),j) < sorted(t)));
        err2 = sum(z(y~=yhat(2)).*(X(y~=yhat(2),j) >= sorted(t)));
        
        err = err1 + err2;
        % If the current split gives less errors, remember it
        if err < minErr && (yhat(1) ~= yhat(2))
            minErr = err;
            minj = j;
            minThreshold = (sorted(t)+sorted(t-1))/2;
            
            % Store predicted class below and above the threshold
            % which class is the most likely in each interval
            predictedClass(1) = yhat(1);
            predictedClass(2) = yhat(2);
        end
    end
end

% Model outputs
model.name = 'Stump Classification';
model.nClasses = nClasses;
model.j = minj;
model.threshold = minThreshold;
model.predictedClass = predictedClass;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);
    
% Split data using threshold
index = Xhat(:,model.j) < model.threshold;

yhat = zeros(nTest, 1);
yhat(index) = model.predictedClass(1);
yhat(~index) = model.predictedClass(2);
end
