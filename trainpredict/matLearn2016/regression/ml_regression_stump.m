function [model] = ml_regression_stump(X,y,options)
% ml_regression_stump(X,y,options)
%
% Description:
%    - Finds the best way to split the training set in two and fits model on
%   each side
%
% Options:
%	 - weights: vector of training data weights (default: vector of 1's)
%	 - error: method used to calculate the minimum error defined by one of:
%       - 'mn' for the absolute/L1 distance from the mean
%       - 'sq' for the squared/L2 distance from the mean
%       (default: 'sq')
%	 - modelType: model used to fit the two subsets of splitted data defined
%       by one of:
%       - 'cns' fitting two constant models using mean of each subset
%       - 'lin' fitting two linear models minizing L2 errors of each subset
%       {default: 'lin')
%    - LambdaL2: strength of L2 regularizaion for linear model (default: 0)
%
% Added Values:
%	 - Select model to calculate error
%	 - Add the possibility to use weights
%
% Authors: Giovanni Viviani (2014)

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
[z, error, lambdaL2, modelType] = myProcessOptions(options,'weights',ones(nTrain, 1),'error','sq', 'lambdaL2', 0, 'modelType', 'lin');

% Constant model
if strcmp(modelType, 'cns');
    
    % Select best feature and splitting threshold to minimize error
    minErr = inf;
    for j = 1:nFeatures
        sorted = sort(unique(X(:,j)));
        for t = 2:length(sorted)
            % Divide training data in two
            ind1 = find(X(:,j) < sorted(t));
            ind2 = find(X(:,j) >= sorted(t));
            % Calculate mean (constant model) on both sides
            mean1 = mean(y(ind1));
            mean2 = mean(y(ind2));
            % Calculate error of each side (accounting for weights)
            switch error
                case 'mn'
                    err1 = z(ind1)'*(abs(y(ind1) - mean1));
                    err2 = z(ind2)'*(abs(y(ind2) - mean2));
                case 'sq'
                    err1 = z(ind1)'*((y(ind1) - mean1).^2);
                    err2 = z(ind2)'*((y(ind2) - mean2).^2);
                otherwise
                    err1 = z(ind1)'*((y(ind1) - mean1).^2);
                    err2 = z(ind2)'*((y(ind2) - mean2).^2);
            end
            % Sum errors from the two sides
            err = err1+err2;
            % Select best error, feature, threshold and constant models
            if err < minErr;
                minErr = err;
                minThreshold = (sorted(t)+sorted(t-1))/2;
                minj = j;
                minmean1 = mean1;
                minmean2 = mean2;
            end
        end
    end
    % Model outputs specific to constant model
    model.modelType = 'cns';
    model.mean1 = minmean1;
    model.mean2 = minmean2;
    
% Linear Model
elseif strcmp(modelType, 'lin');
    % Add bias variable
    XwithBias = [ones(nTrain, 1), X];
    % Select best feature and splitting threshold to minimize error
    minErr = inf;
    for j = 1:nFeatures
        thresholds = [sort(unique(X(:,j)));max(X(:,j))+eps]';
        for t = thresholds
            % Divide training data in two
            ind1 = find(X(:,j) < t);
            ind2 = find(X(:,j) >= t);
            % Calculate linear model on both sides
            [w1, yhat1] = leastSquareL2normWeight(XwithBias(ind1, :), y(ind1), z(ind1), lambdaL2);
            [w2, yhat2] = leastSquareL2normWeight(XwithBias(ind2, :), y(ind2), z(ind2), lambdaL2);
            % Calculate error of each side (accounting for weights)
            switch error
                case 'mn'
                    err1 = z(ind1)'*(abs(y(ind1) - yhat1));
                    err2 = z(ind2)'*(abs(y(ind2) - yhat2));
                case 'sq'
                    err1 = z(ind1)'*((y(ind1) - yhat1).^2);
                    err2 = z(ind2)'*((y(ind2) - yhat2).^2);
                otherwise
                    err1 = z(ind1)'*((y(ind1) - yhat1).^2);
                    err2 = z(ind2)'*((y(ind2) - yhat2).^2);
            end
            % Sum errors from the two sides
            err = err1+err2;
            % Select best error, feature, threshold and linear models
            if err < minErr;
                minErr = err;
                minThreshold = t;
                minw1 = w1;
                minw2 = w2;
                minj = j;
            end
        end
    end
    % Model outputs specific to linear model
    model.modelType = 'lin';
    model.w1 = minw1;
    model.w2 = minw2;
end
% Model outputs for both models
model.name = 'Stump Regression';
model.threshold = minThreshold;
model.predict = @predict;
model.j = minj;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
yhat = zeros(nTest,1);

index1 = Xhat(:,model.j) < model.threshold;
index2 = Xhat(:,model.j) >= model.threshold;
    
% Constant Model
if strcmp(model.modelType, 'cns');
    % Select usage between the two submodels for each test data
    yhat(index1) = model.mean1;
    yhat(index2) = model.mean2;
end
% Linear Model
if strcmp(model.modelType, 'lin');
    % Add bias variable
    X = [ones(nTest,1) Xhat];
    % Select usage between the two submodels for each test data
    yhat(index1) = X(index1,:)*model.w1;
    yhat(index2) = X(index2,:)*model.w2;
end
end

function [w yhat] = leastSquareL2normWeight(X, y, z, lambda)
% Creates linear regression (L2) model with training sample weights and L2
%   regularization

% Optimization options for minFunc
optimOptions.Display = 0; 
optimOptions.useMex = 0;

% compute weight vector w using minFunc
w = minFunc(@squaredLossL2,randn(size(X,2),1),optimOptions,X,y,lambda,z);
yhat = X * w;
end

function [f,g] = squaredLossL2(w,X,y,lambda,z)
% Squared loss function with training sample weights and L2 regularization
f = (X*w-y)'*diag(z)*(X*w-y) + lambda*(w'*w);
g = 2*X'*diag(z)*(X*w-y) + 2*lambda*w;
end
