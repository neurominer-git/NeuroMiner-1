function [model] = ml_generative_NB(X,~,options)
% ml_multiclass_generativeNB(X,~,options)
%
% Description:
%    - Computes the maximum likelihood Naive Bayes (Gaussian/categorical counts)
%       fit for the data.Targets y are expected to be empty and are ignored.
%       The prediction function returns an [nTest, 1] vector with the
%       likelihood of each data point in the test set.
%
% Options:
%    - dataTypes: vector 0/1's where 0 indicate i'th feature is discrete,
%       and 1 means feature is continuous (default: vector of 1's)
%
% Authors:
% 	 - Celia Siu (2014)

[nTrain, nFeatures] = size(X);

% Default options
[dataTypes] = myProcessOptions(options,'dataTypes', ones(nFeatures,1));

% Calculate and save 1-dimensional distribution data for each feature
for j = 1:nFeatures
    if ~dataTypes(j)
        % Fraction of counts of data in each category
        cats = unique(X(:,j));
        numCat = size(cats, 1);
        for k = 1:numCat
            fcounts = sum(X(:,j) == cats(k));
            features{j}.prob(k) = fcounts / nTrain;
        end
    elseif dataTypes(j)
        % 1-dimensional Gaussian parameters
        features{j}.mu = mean(X(:,j));
        features{j}.sigma = sqrt(var(X(:,j)));
    end
end

% Model outputs
model.name = 'Generative Naive Bayes Model';
model.dataTypes = dataTypes;
model.features = features;
model.predict = @predict;

end

function [lik] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

% Multiply distributions of each feature
lik = 1;
for j = 1:nFeatures
    if ~model.dataTypes(j)
        try
            lik = lik .* model.features{j}.prob(Xhat(:,j));
        catch
            % When test set feature label is not in training set
            lik = lik * 0;
        end
    elseif model.dataTypes(j)
        lik = lik .* normpdf(Xhat(:,j),model.features{j}.mu, model.features{j}.sigma);
    end
end
end
