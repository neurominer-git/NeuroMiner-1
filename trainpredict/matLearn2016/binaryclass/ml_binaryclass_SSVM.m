function [model] = ml_binaryclass_SSVM(X,y,options)
% ml_binaryclass_SSVM(X,y,options)
%
% Description:
%	 - Fits a smooth support vector machine by minimizing squared hinge loss
%
% Options:
%    - addBias: adds a bias variable (default: 0)
%    - lambdaL2: strength of L2-regularization parameter (default: 0)
%    - kernel: use kernel function (default: 0)
%    - kernelFunc: kernel function (default: rbf)
%    - kernelOptions: options passed to kernel function
%       (default: sigma = 1)
%
% Note: The alogrithms below were based on the examples posted on Professor
%       Schmidt's website. Some functions were taken directly from here
%       http://www.cs.ubc.ca/~schmidtm/Software/minFunc/examples.html#16
%
% Authors:
%  	 - Manyou Ma (2014)

[nTrain,nFeatures] = size(X);

% Use input options if they exist, otherwise use defaults
kernelFuncOptions = [];
kernelFuncOptions.sigma = 1;
[addBias,kernel,kernelFunc,kernelOptions,lambdaL2] = myProcessOptions(options,'addBias',0,'kernel',0, 'kernelFunc', @ml_kernel_rbf, 'kernelOptions', kernelFuncOptions, 'lambdaL2', 1.0);

% Add the bias variable
if (addBias==1)
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Apply kernelization
if kernel
    model.base = X;
    X = kernelFunc(X,X,kernelOptions);
    nFeatures = nTrain;
    model.newX  = X;
    
end

% Optimization options
optimOptions.Display = 0;
optimOptions.useMex = 1;

% Find the optimal estimation using optimization
if lambdaL2 == 0
    % Without regularization
    funObj = @(w)SSVMLoss(w,X,y);
    w = minFunc(funObj,zeros(nFeatures,1),optimOptions);
else
    % L2 regularization
    funObjL2 = @(w)SSVMLossL2(w,X,y,lambdaL2);
    w = minFunc(funObjL2,zeros(nFeatures,1),optimOptions);
end

% Model output
model.name = 'Squared Hinge Loss SVM Binary Classification';
model.kernel = kernel;
model.kernelFunc = kernelFunc;
model.kernelOptions = kernelOptions;
model.addBias = addBias;
model.w = w;
model.predict = @predict;
end

function p = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);

% Optional bias
if (model.addBias==1)
    Xhat = [ones(nTest,1) Xhat];
end

% Optional kernalize
if (model.kernel == 1)
    newX = model.kernelFunc(Xhat,model.base,model.kernelOptions);
    Xhat = newX;
end

p.D = Xhat*model.w; p.D(isnan(p.D))=0; 
index =  p.D >= 0;
p.yhat = sign(index - 0.5);

end

function [f,g] = SSVMLoss(w,X,y)
% Smooth squared hinge loss function
% Taken from Prof. Schmidt's website
err = 1-y.*(X*w);
viol = find(err>=0);
f = sum(err(viol).^2);

if isempty(viol)
    g = zeros(size(w));
else
    g = -2*X(viol,:)'*(err(viol).*y(viol));
end
end

function [f,g] = SSVMLossL2(w,X,y,lambda)
% Smooth squared hinge loss function with L2 regularization
err = 1-y.*(X*w);
viol = find(err>=0);
f = sum(err(viol).^2);
f = f+sum(lambda.*(w.^2));

if isempty(viol)
    g = zeros(size(w));
else
    g = -2*X(viol,:)'*(err(viol).*y(viol));
end
g = g + 2*lambda.*w;
end