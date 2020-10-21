function [model] = ml_kernel_ordinal_logistic(X,y,options)
% [model] = ml_kernel_ordinal_logistic(X,y,options)
%
% Description:
%    - Classification using kernel ordinal logistic regression model
%
% Options:
%   - verbose: logical 0/1 flag controlling verbosity of optimization
%               algorithm (default: 1)
%   - lambdaL2: scalar constant regularizer for all weights (default: 1e-5)
%   - kernelFunc: kernel function handle for ordinal logistic regression
%                 (default: @ml_kernel_poly) 
%   - kernelOptions: struct containing options for kernelFunc (default:
%                    empty)
%
% Authors:
% 	- Mark Schmidt (2014)

[n,p] = size(X);

if nargin < 3
    options = [];
end

[verbose,lambdaL2,kernelFunc,kernelOptions] = myProcessOptions(options, ...
   'verbose',0,'lambdaL2',1e-5,'kernelFunc',@ml_kernel_poly, ...
   'kernelOptions', []);
k = options.nClasses;

if verbose
    optimoptions.verbose = 3;
else
    optimoptions.verbose = 0;
end

% Form nxn Gram matrix K
[K, kernName] = kernelFunc(X,X,kernelOptions);

% Initialize thresholds to learn
u = zeros(n,1);
gamma = ones(k-2,1);
LB = [-inf(n,1);zeros(k-2,1)];
UB = inf(n+k-2,1);

funObj_sub = @(u)ml_ordinal_logistic_loss(u,K,y,k);
% Select which weights to be L2-regularized
subset = 1:n;
funObj = @(u)ml_penalized_kernel_L2_subset(u,K,subset,funObj_sub,lambdaL2);

% Do weights and threshold optimization
u_gamma = minConf_TMP(funObj,[u(:);gamma(:)],LB,UB,optimoptions);
u = u_gamma(1:n);
gamma = [-inf;0;cumsum(u_gamma(n+1:end));inf];

model.name = strcat(['Ordinal Logistic Regression, ',kernName]); 
model.nClasses = k;
model.thresh = gamma;
model.predict = @predict;
model.lossFunc = @loss;
model.weights = u;
model.Xtrain = X;
model.kernelFunc = kernelFunc;
model.kernelOptions = options.kernelOptions;
end

function y = predict(model,X)
k = model.nClasses;
u = model.weights;
gamma = model.thresh;
[n,p] = size(X);
z = model.kernelFunc(X,model.Xtrain,model.kernelOptions)*u;
y = zeros(n,1);
for c = 1:k
    y(z > gamma(c)) = c;
end
end

function f = loss(model,X,y)
u = model.weights;
gamma = model.thresh;
K = model.kernelFunc(X,model.Xtrain,model.kernelOptions);
z = K*u;
sigmoid1 = 1./(1+exp(z - gamma(y+1)));
sigmoid2 = 1./(1+exp(z - gamma(y)));
f = -sum(log(sigmoid1-sigmoid2+eps))';
end
