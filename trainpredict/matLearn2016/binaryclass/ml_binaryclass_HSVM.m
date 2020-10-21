function [model] = ml_binaryclass_HSVM(X,y,options)
% ml_binaryclass_HSVM(X,y,options)
%
% Description:
%   - Fits a linear classifier using Huberized SVM
%
% Options:
%	 - addBias: adds a bias variable (default: 0)
%	 - epsilon: threshold for Huberized hinge loss (default: 0.5)
%    - lambdaL2: strength of L2-regularization parameter (default: 0)
%    - lambdaL1: strength of L1-regularization parameter (default: 0)
%       * if both lambdaL2 & lambdaL1 are provided, use Elastic Net Penalty
%    - kernel: apply kernel (default: 0)
%    - kernelFunc: kernel function used (default: rbf)
%    - kernelOptions: options passed to kernel function
%       (default: sigma = 1)
%
% Notes:
%   - The L1-regularization uses the solver from package L1General
%   - This function is completed with the help from L1GeneralExamples on
%   Professor Schmidt's website, http://www.cs.ubc.ca/~schmidtm/Software/L1General.html
%
% Authors:
% 	- Ben Zhu (2014)

[nTrain,nFeatures] = size(X);

% Process options
kernelFuncOptions = [];
kernelFuncOptions.sigma = 1;
[addBias,epsilon,lambdaL2,lambdaL1, kernel, kernelFunc, kernelOptions] = ...
    myProcessOptions(options,'addBias',0,'epsilon',0.5,'lambdaL2',0,'lambdaL1',0,'kernel',0,'kernelFunc', @ml_kernel_rbf, 'kernelOptions', kernelFuncOptions);

% Add the bias variable
if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Kernalize if required
if kernel
    model.base = X;
    X = kernelFunc(X,X,kernelOptions);
    nFeatures = nTrain;
    model.newX  = X;
end

% Optimization options
optimOptions.Display = 0;
optimOptions.useMex = 1;

% Fit in solvers
if lambdaL2 == 0 && lambdaL1 == 0
    % No regularization
    funObj = @(w)HSVMLoss(w,X,y,epsilon);
    w = minFunc(funObj,zeros(nFeatures,1),optimOptions);
elseif lambdaL1 == 0
    % L2 regularization
    funObjL2 = @(w)HSVMLossL2(w,X,y,epsilon,lambdaL2);
    w = minFunc(funObjL2,zeros(nFeatures,1),optimOptions);
elseif lambdaL2 == 0
    % L1 regularization
    funObj = @(w)HSVMLoss(w,X,y,epsilon);
    lambda = lambdaL1*ones(nFeatures,1);
    % Do not penalize bias variable
    lambda(1) = 0;
    w = L1General2_PSSgb(funObj,zeros(nFeatures,1),lambda);
else
    % Elastic net penalty
    lambda = lambdaL1*ones(nFeatures,1);
    % Do not penalize bias variable
    lambda(1) = 0;
    funObjL2 = @(w)HSVMLossL2(w,X,y,epsilon,lambdaL2);
    w = L1General2_PSSgb(funObjL2,zeros(nFeatures,1),lambda);
end
if numel(w) ~= nFeatures
    fprintf('W problem');
end
% Model outputs
model.name = ['Huberized Hinge SVM Binary Classification, ', num2str(epsilon), ' \epsilon'];
model.addBias = addBias;
model.kernel = kernel;
model.kernelFunc = kernelFunc;
model.kernelOptions = kernelOptions;
model.w = w;
model.predict = @predict;
end

function p = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);

% Add optional bias
if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end

% Kernalize
if model.kernel
    newX = model.kernelFunc(Xhat,model.base,model.kernelOptions);
    Xhat = newX;
end

p.D = Xhat*model.w; p.D(isnan(p.D))=0; 
index =  p.D >= 0;
p.yhat = sign(index - 0.5);
end

function [f,g] = HSVMLoss(w,X,y,epsilon)
% Huberized hinge loss from L1GeneralExamples
[nTrain,nFeatures] = size(X);
f = 0;
g = zeros(nFeatures,1);
yhat = y.*(X*w);
ind1 = yhat <= epsilon;
ind2 = yhat > epsilon & yhat <= 1;
if any(ind1)
    f = f + sum((1-epsilon)^2 + 2*(1-epsilon)*(epsilon-yhat(ind1)));
    g = g - X(ind1,:)'*(2)*(1-epsilon)*y(ind1);
end
if any(ind2)
    f = f + sum((1-yhat(ind2)).^2);
    g = g - X(ind2,:)'*2*((1-yhat(ind2)).*y(ind2));
end
end

function [f,g] = HSVMLossL2(w,X,y,epsilon,lambda)
% Huberized hinge loss with L2 regularization
[nTrain,nFeatures] = size(X);
f = 0;
g = zeros(nFeatures,1);
yhat = y.*(X*w);
ind1 = yhat <= epsilon;
ind2 = yhat > epsilon & yhat <= 1;
if any(ind1)
    f = f + sum((1-epsilon)^2 + 2*(1-epsilon)*(epsilon-yhat(ind1)));
    g = g - X(ind1,:)'*(2)*(1-epsilon)*y(ind1);
end
if any(ind2)
    f = f + sum((1-yhat(ind2)).^2);
    g = g - X(ind2,:)'*2*((1-yhat(ind2)).*y(ind2));
end
f = f + (lambda/2)*(w'*w);
g = g + lambda*w;
end