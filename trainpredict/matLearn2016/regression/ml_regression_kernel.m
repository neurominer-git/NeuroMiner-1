function [model] = ml_regression_kernel(X,y,options)
% ml_regression_kernel(X,y,options)
%
% Description:
%	 - Fitting a linear regression model by minimizing the squared loss of a kernel
%	   
% Options:
%	 - weights: vector of training data weights (default: vector of 1's)
%	 - lambdaL2: strength of L2-regularization parameter (default: 1)
%	 - addBias: accepts 0 or 1. If 1, adds bias to X (default: 1)
%   	 - kernelFunc: kernel function for kernel trick (default: Gram Matrix)
%   	 - kernelOptions: options for kernel function (default: none)

[nTrain,nFeatures] = size(X);

% Generate default model options if not given by input
[lambdaL2,addBias,kernelFunc,nClasses,kernelOptions] = myProcessOptions(options,...
    'lambdaL2',1,'addBias',1,'kernelFunc',@ml_kernel_gram,'nClasses',...
    0,'kernelOptions', []);

% Must have L2 regularization
if lambdaL2 == 0
    fprintf('ERROR: L2 regularization must be applied.');
    return;
end

% Add bias variable if required
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

[kern,name] = kernelFunc(X,X, kernelOptions);
a = (kern + lambdaL2 *eye(nTrain))\y;

% Model outputs
model.name = ['LS Regression, ', name];
model.nClasses = nClasses;
model.a = a;
model.X = X;
model.addBias = addBias;
model.kernelFunc = kernelFunc;
model.kernelOptions = kernelOptions;
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
% Add optional bias variable
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
    nFeatures = nFeatures + 1;
end
yhat = model.kernelFunc(model.X, Xhat, model.kernelOptions)'*model.a;
end
