function [model] = ml_regression_GAM(X,y,options)
% ml_regression_GAM(X,y,options)
%
% Description:
%	 - Fits a general additive model using backfitting algorithm
%
% Options:
%    - subFunc: function used to determine f's (default: smooth cubic splines)
%       - 'spl' - smooth cubic splines
%       - 'rg' - polynomial regression
%       - 'lin' - linear regresion
%    - deg: highest exponent degree for polynomial subFunc (default: 2)
%    - smoothing: smoothness of smooth cubic splines subFunc, where 0 is
%       linear regression and 1 is unsmoothed cublic splines (default: 0.8)
%    - maxIter: max number of iteratons for backfitting (default: 500)
%
% Author: Jennifer (Xin Bei) She (2015)

[nTrain,nFeatures] = size(X);

% Default values
[maxIter, subFunc, deg, smoothing, maxIterS] = ...
    myProcessOptions(options,'numIter',500, 'subFunc', 'spl', 'deg', 2, ...
                     'smoothing', 0.8, 'maxIter', 500);

% Backfitting algorithm
% Initialize alpha
alpha = (1/nTrain)*sum(y);
offset = zeros(nFeatures, 1);

% Initialize f
for i=1:nFeatures
    if strcmp(subFunc, 'sp1')
        f{i} = csaps(X(:,i), zeros(nTrain, 1), smoothing);
    elseif strcmp(subFunc, 'rg')
        f{i} = polyfit(X(:,i), zeros(nTrain, 1), deg);
    elseif strcmp(subFunc, 'lin')
        f{i} = ml_regression_L2(X(:,i), zeros(nTrain, 1), []);
    end
end

prevError = inf;
flag = 0;
% Cycle until error starts increasing (unless # iterations exceeds maxIter)
for i=1:maxIter
    for j=1:nFeatures
        
        % Calculate yhat
        yhat = alpha*ones(nTrain, 1);
        for k=1:nFeatures
            if strcmp(subFunc, 'spl')
                yhat = yhat + fnval(f{k}, X(:,k)) - offset(k);
            elseif strcmp(subFunc, 'rg')
                yhat = yhat + polyval(f{k}, X(:,k)) - offset(k);
            elseif strcmp(subFunc, 'lin')
                yhat = yhat + f{k}.predict(f{k}, X(:,k)) - offset(k);
            end
            
        end
        
        % Break if error starts increasing
        currError = sum((y-yhat).^2);
        if prevError <= currError
            flag = 1;
            break;
        else
            prevError = currError;
        end
        
        % Update fj
        if strcmp(subFunc, 'spl')
            yhat = yhat - fnval(f{j}, X(:,j)) + offset(j);
            f{j} = csaps(X(:,j), y-yhat, smoothing);
            offset(j) = (1/nTrain)*sum(fnval(f{j}, X(:, j)));
        elseif strcmp(subFunc, 'rg')
            yhat = yhat - polyval(f{j}, X(:,j)) + offset(j);
            f{j} = polyfit(X(:,j), y-yhat, deg);
            offset(j) = (1/nTrain)*sum(polyval(f{j}, X(:, j)));
        elseif strcmp(subFunc, 'lin')
            yhat = yhat - f{j}.predict(f{j}, X(:,j)) + offset(j);
            f{j} = ml_regression_L2(X(:,j), y-yhat, []);
            offset(j) = (1/nTrain)*sum(f{j}.predict(f{j}, X(:, j)));
        end
    end
    if flag
        break;
    end
end

% Model outputs
model.f = f;
model.XTrain = X;
model.nTrain = nTrain;
model.predict = @predict;
model.alpha = alpha;
model.offset = offset;
model.subFunc = subFunc;
if strcmp(subFunc, 'spl')
    model.name = 'GAM with Smooth Cubic Splines';
elseif strcmp(subFunc, 'rg')
    model.name = ['GAM with Polynomial Regression of Degree ', num2str(deg)];
elseif strcmp(subFunc, 'lin')
    model.name = 'GAM with Linear Regression';    
end
end

function [yhat] = predict(model, Xhat)
% Predict function
[nTest,nFeatures] = size(Xhat);

% Add up respective f(Xhat)'s
yhat = model.alpha*ones(nTest, 1);
if strcmp(model.subFunc, 'spl')
    for i=1:nFeatures
        yhat = yhat + fnval(model.f{i}, Xhat(:,i)) - model.offset(i);
    end
elseif strcmp(model.subFunc, 'rg')
    for i=1:nFeatures
        yhat = yhat + polyval(model.f{i}, Xhat(:,i)) - model.offset(i);
    end
elseif strcmp(model.subFunc, 'lin')
    for i=1:nFeatures
        yhat = yhat + model.f{i}.predict(model.f{i}, Xhat(:,i)) - model.offset(i);
    end
    
end
end