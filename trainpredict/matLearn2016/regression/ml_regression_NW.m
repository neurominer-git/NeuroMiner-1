function [model] = ml_regression_NW(X,y, options)
% ml_regression_NW(X,y, options)
%
% Description:
%    - Uses a kernel function to do a weighted average KNN-style regression
%
% Options:
%    - kernelFunc: a kernel function that calculates a weight given the
%       difference between two vectors (default: Gaussian)
%    - kernelOptions: options for kernel fuction
%       (default: kernel variance = 1)
%
% Authors:
%    - Ben Bougher (2014)

% Store the data
kernel_options = [];
kernel_options.sigma = 0.5;
[kernelFunc, kernelOptions] = myProcessOptions(options,'kernelFunc', @ml_kernel_rbf, 'kernelOptions', kernel_options);

% Store the model outputs
model.X = X;
model.y = y;
[~, name] =kernelFunc(X,X,options);
model.name = ['K-Nearest Neighbour Regression Weighted with: ', name];
model.predict = @predict;
model.kernelFunc = kernelFunc;
model.kernelOptions =kernelOptions;
end

function [yhat] = predict(model, Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
[nTrain, nFeatures] = size(model.X);
yhat = zeros(nTest,1);

% Loop through the prediction and training data
for i = 1:nTest
    n = zeros(nTrain, 1);
    % Calculate amount of contribution of each neighbour
    for j = 1:nTrain
        n(j) = model.kernelFunc(Xhat(i,:),model.X(j,:),model.kernelOptions);
    end
    % Caculate weighted average 
    yhat(i) = sum(model.y.*n)/sum(n);
end
end
