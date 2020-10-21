function [model] = ml_regression_SVR(X,y,options)
% ml_regression_SVR(X,y,options)
%
% Description:
%	 - Support vector regression using the eps-insensitive loss
%
% Options:
%    - addBias: add a bias variable for un-kernalized SVR (default: 1)
%    - epsilon: zero error if the absolate difference between
%       the prediction wx_i and the target y_i is less than eps
%       (default: 1, problem becomes a ridge regression)
%    - lambdaL2: regularization coefficient (default: 0)
%    - method: implementation of minimizing error, one of:
%       'sgd' - solve the primal problem with SGD
%       'sm' - solve the smoothed primal quadratic
%       programming problem
%       'smk' - solve the smoothed primal quadratic
%        programming problem, kernelized
%       (default: 'sm')
%    - kernelFunc: kernel function (default: polyKernel)
%    - kernelOptions: (default: bias = 1, order = 1)
%    - threshold: use |w_t+1 - w_t| < threshold as the stopping criteria
%       (default: 1e-5)
%    - maxIter: use maximum iterations as the stopping criteria
%       (default: 10000)
%    * when both threshold and maxIter are provided, threshold will be
%       pursued, until maximum iteration is hit
%
% Authors:
% 	 - Yan Peng (2014)
%

% Setting up parameters
kernelDefaultOptions = [];
kernelDefaultOptions.sigma = 1;
[method,kernelFunc,kernelOptions,lambdaL2,epsilon,maxIter,threshold,...
addBias] = myProcessOptions(options,'method','sgd',...
                            'kernelFunc',@ml_kernel_rbf,...
                            'kernelOptions', kernelDefaultOptions,...
                            'lambdaL2',0,...
                            'epsilon',1.0,...
                            'maxIter',10000,...
                            'threshold',1e-5,...
                            'addBias',1);

[nTrain, nFeatures] = size(X);

% Add bias variable if required
if options.addBias || strcmp(method, 'smk');
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Setup parameters
options.method = method;
options.lambda = lambdaL2;
options.epsilon = epsilon;
options.addBias = addBias;

if (strcmp(method,'sgd'))
    % Set up maxIter and threshold for SGD
    options.maxIter = maxIter;
    options.threshold = threshold;
    % Train primal problem with SGD
    [w] = SGDPrimal(X,y, options);
    model.supportVector = abs(X*w - y) >= epsilon;
elseif (strcmp(method,'sm'))
    % Train primal problem with smoothness
    [w] = SmoothPrimal(X,y, options);
    model.supportVector = abs(X*w - y) >= epsilon;
elseif (strcmp(method,'smk'))
    % Set up kernel function
    options.kernelFunc = kernelFunc;
    options.kernelOptions = kernelOptions;
    % Train primal problem with smoothness, kernelized
    [w, K] = SmoothPrimalKernel(X,y, options);
    model.kernelFunc = kernelFunc;
    model.kernelOptions = kernelOptions;
    model.supportVector = abs(K*w - y) >= epsilon;
    model.X2 = X(:,2);
end

% Model outputs
model.name = ['Support Vector Regression with ',num2str(epsilon),' Epsilon'];
model.w = w;
model.epsilon = options.epsilon;
model.predict = @predict;
model.addBias = options.addBias;
model.method = method;
model.kernelFunc = kernelFunc;
end

function [w] = SGDPrimal(X,y, options)
% Stochastic gradient descent minimization
[nTrain, nFeatures] = size(X);

% Initial value of w
w = zeros(nFeatures,1);

% Evaluate function and gradient at initial point
f = sum(max(0,abs(y-X*w)-options.epsilon)) + options.lambda*(w'*w)/2;

% Matlab stores sparse vectors as columns, so when accessing individual
% rows it is much faster to work with X^T
Xt = X';

% Optimize primal problem using stochastic gradient descent method
iter = 0;
w_old = ones(nFeatures,1)*Inf;
% Check if the difference between successive w is small enough
% Also iter < p.maxIter/10 to ensure at least that many iterations
while abs(sum(w_old - w)) > options.threshold || iter < options.maxIter/10
    iter = iter + 1;
    w_old = w;
    
    % Choose a random integer between 1 and N
    i = ceil(nTrain*rand);
    
    % Compute a subgradient with respect to example i
    if (abs((Xt(:,i)'*w)-y(i)) <= options.epsilon)
        g = options.lambda*w;
    elseif ( y(i)-(Xt(:,i)'*w) > 0 )
        g = -Xt(:,i)+options.lambda*w;
    else
        g = +Xt(:,i)+options.lambda*w;
    end
    
    % Step size & params update
    alpha = 1/iter;
    w = w - alpha*g;
    
    % Check if exceed maximum iterations
    if (iter > options.maxIter)
        fprintf('Maximum iterations...\n');
        break;
    end
end

% Evaluate function and gradient at ending point
f = sum(max(0,abs(y-X*w)-options.epsilon)) + options.lambda*(w'*w)/2;
end

function [w, K] = SmoothPrimal(X,y,options)
% Minimization using quadratic programming
[nTrain, nFeatures] = size(X);

% Setting up linear/quadratic programming parameters
f = [zeros(nFeatures,1);ones(nTrain,1)];
A = [-X, -eye(nTrain); ...
     X, -eye(nTrain); ...
     zeros(nTrain,nFeatures), -eye(nTrain)];
b = [-y+options.epsilon;y+options.epsilon;zeros(nTrain,1)];

if options.lambda == 0
    % Linear programming when lambda is 0
    options_ = optimoptions(@linprog,'MaxIter',5000,...
        'Display', 'off');
    [wv] = linprog(f,A,b,[],[],[],[],[],options_);
else
    % Quadratic programming otherwise
    options_ = optimoptions(@quadprog,'MaxIter',5000,...
        'Display', 'off');
    H = [2*options.lambda*eye(nFeatures),zeros(nFeatures,nTrain);...
        zeros(nTrain,nFeatures),zeros(nTrain,nTrain)];
    [wv] = quadprog(H,f,A,b,[],[],[],[],[],options_);
end
w = wv(1:nFeatures);
end

function [w, K] = SmoothPrimalKernel(X,y,options)
% Minimization using quadratic programming, kernalized
[nTrain, nFeatures] = size(X);

% Use kernel matrix instead
K = options.kernelFunc(X,X,options.kernelOptions);

% Setting up linear/quadratic programming parameters
f = [zeros(nTrain,1);ones(nTrain,1)];
A = [-K, -eye(nTrain); ...
    K, -eye(nTrain); ...
    zeros(nTrain,nTrain), -eye(nTrain)];
b = [-y+options.epsilon;y+options.epsilon;zeros(nTrain,1)];

if options.lambda == 0
    % Linear programming when lambda is 0
    options_ = optimoptions(@linprog,'MaxIter',5000,...
        'Display', 'off');
    [wv] = linprog(f,A,b,[],[],[],[],[],options_);
else
    % Quadratic programming otherwise
    options_ = optimoptions(@quadprog,'MaxIter',5000,...
        'Display', 'off');
    H = [2*options.lambda*eye(nFeatures),zeros(nFeatures,nTrain);...
        zeros(nTrain,nFeatures),zeros(nTrain,nTrain)];
    [wv] = quadprog(H,f,A,b,[],[],[],[],[],options_);
end
w = wv(1:nTrain);
end

function [yhat] = predict(model,Xhat)
% Prediction variable
[nTest, nFeatures] = size(Xhat);

% Kernalized
if strcmp(model.method, 'smk')
    yhat = model.kernelFunc(Xhat, model.X2, model.kernelOptions)*model.w;
else
    % Un-kernealized
    % Add a bias variable
    if model.addBias;
        Xhat = [ones(nTest,1) Xhat];
        nFeatures = nFeatures + 1;
    end
    yhat = Xhat*model.w;
end
end