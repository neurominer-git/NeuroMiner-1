function [model] = ml_multiclass_kernel_softmax(X,y,options)
% [model] = ml_classification_kernel_softmax(X,y,options)
%
% Description:
%	 - Fits a classification model by minimizing the kernel softmax 
%      loss function
%
% Options:
%   - verbose: logical 0/1 to toggle verbosity level 
%   - kernelFunc: kernel function (default: Linear kernel)
%   - kernelOptions: options for kernel function (default: none)
%	- lambdaL2: strength of L2-regularization parameter (default: 1)
%	- link: 'ssvm' for structured SVM or 'softmax' for softmax loss



[nInstances,nFeatures] = size(X);

if nargin < 3
    options = [];
end

[verbose,kernelFunc,kernelOptions,lambdaL2,link] = ...
    myProcessOptions(options,'verbose',0,'kernelFunc',@ml_kernel_gram, ...
    'kernelOptions',options.kernelOptions,'lambdaL2',1e-5,'link','softmax');

if verbose
    optimoptions.Display = 'iter';
else
    optimoptions.Display = 'none';
end

k = options.nClasses;
[K, kernName] = kernelFunc(X,X,kernelOptions);
gradArgs = {K,y,k};

% L2-Regularization
if strcmp(link,'ssvm')
    lossFunc = @ml_SSVM_multiLoss;
    u = zeros(nInstances,k);
    u(:) = minFunc(@ml_penalized_kernel_L2_matrix,u(:),optimoptions,K,k, ...
                   lossFunc,lambdaL2,gradArgs{:});
else
    lossFunc = @ml_softmax_loss;
    u = zeros(nInstances,k-1);
    u(:) = minFunc(@ml_penalized_kernel_L2_matrix,u(:),optimoptions, ...
                   K,k-1,lossFunc,lambdaL2,gradArgs{:});
end
model.name = strcat([kernName,' ','Classification']);
model.nTrain = nInstances;
model.nClasses = k;
model.weights = u;
model.Xtrain = X;
model.kernelFunc = kernelFunc;
model.kernelOptions = options.kernelOptions;
if strcmp(link,'ssvm')
    model.name = strcat([model.name,', SSVM Loss']);
    model.predict = @predictSVM;
else
    model.name = strcat([model.name,', Softmax Loss']);
    model.predict = @predictSoftmax;
end
model.lossFunc = @(model,X,y)ml_softmax_loss(model.weights,...
                                    model.kernelFunc(X, model.Xtrain,...
                                                     options.kernelOptions),...
                                    y,model.nClasses);
end

function y = predictSoftmax(model,X)
k = model.nClasses;
u = model.weights;
X = model.kernelFunc(X,model.Xtrain,model.kernelOptions);
[n,p] = size(X);
[junk y] = max(X*[u zeros(model.nTrain,1)],[],2);
end

function y = predictSVM(model,X)
k = model.nClasses;
u = model.weights;
X = model.kernelFunc(X,model.Xtrain,model.kernelOptions);
[n,p] = size(X);
[junk y] = max(X*u,[],2);
end
