function [model] = ml_binaryclass_MLP(X,y,options)
% ml_binaryclass_MLP(X,y,options)
%
% Description:
%	 - Binary classification using a multilayer perceptron with logistic
%       loss
%
% Options:
%    - nHidden: a vector of the size of each hidden layer,
%       including bias variable (default: [3,3,3])
%    - standardize: if set to 1, standardize columns of X to have
%       mu 0, sigma 1 (default: 0)
%
% Authors:
% 	- Jimmy Yi Pin Cao (2014)

[nTrain,nFeatures] = size(X);

% process inputs
[standardize,nHidden] =...
    myProcessOptions(options,...
    'standardize', 0, 'nHidden', 3*ones(1,3));

% Optimalization options
optimOptions.Display = 0;
optimOptions.useMex = 1;

if numel(nHidden)==1, nHidden = repmat(nHidden,1,nHidden); end

% Standardize optionally
if standardize
    mu = mean(X);
    sigma = std(X);
    X = X - repmat(mu, nTrain, 1);
    X = X ./ repmat(sigma, nTrain, 1);
    % make sure we can apply same transformation at testing time.
    model.standardize = 1;
    model.mu = mu;
    model.sigma = sigma;
end

% Initialize weights
bound = sqrt(6/(nFeatures+nHidden(1)));
w = 2*bound*rand(nFeatures*nHidden(1),1) - bound;

for i = 2:length(nHidden);
    bound = sqrt(6/(nHidden(i-1)+nHidden(i)));
    w = [w; 2*bound*rand(nHidden(i-1)*nHidden(i),1) - bound];
end
bound = sqrt(6/(nHidden(end)+1));
w = [w;  2*bound*rand(nHidden(end),1) - bound];

% Optimalization options
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Optimize weights
w = minFunc(@MLP,w,optimOptions,X,y,nHidden);

% Model outputs
model.nHidden = nHidden;
model.name = 'Multi-layer Perceptron Binary Classification';
model.w = w;
model.predict = @predict;
model.standardize = standardize;
end

function [f,g] = MLP(w, X, y, nHidden)
% Loss function for MLP

% Check for improper input
if(any(isnan(w))); f=inf; g=zeros(size(w)); return; end;

% Get number of instances, variables, and layers
[nTrain,nFeatures] = size(X);

% Extract weights from weight vector
if (~isempty(nHidden))
    inputWeights = reshape(w(1:nFeatures*nHidden(1)),nFeatures,nHidden(1));
    offset = nFeatures*nHidden(1);
    for h = 2:length(nHidden)
        hiddenWeights{h-1} = reshape(w(offset+1:offset+nHidden(h-1)*nHidden(h)),nHidden(h-1),nHidden(h));
        offset = offset+nHidden(h-1)*nHidden(h);
    end
else
    offset = 0;
end
outputWeights = w(offset+1:offset+nHidden(end));

% Initialize function and gradiant
if nargout > 1
    if(~isempty(nHidden))
        gInput = zeros(size(inputWeights));
    else
        gInput = [];
    end
    gOutput = zeros(size(outputWeights));
    for h = 1:length(nHidden)-1
        gHidden{h} = zeros(size(hiddenWeights{h}));
    end
end

% Compute Output
if(~isempty(nHidden))
    ip{1} = X*inputWeights;
    fp{1} = tanh(ip{1});
    
    % Correct bias unit
    ip{1}(:,1) = -inf;
    fp{1}(:,1) = 1;
    for h = 2:length(nHidden)
        ip{h} = fp{h-1}*hiddenWeights{h-1};
        fp{h} = tanh(ip{h});
        
        %Correct for bias unit
        ip{h}(:,1) = -inf;
        fp{h}(:,1) = 1;
    end
    prob = fp{end}*outputWeights;
else
    prob = X*outputWeights{1};
end

% Compute error matrix
f = sum(log(1 + exp(-y.*prob)));

% Compute gradient
if nargout > 1
    
    % Output weight gradient
    err = -(y./(1+exp(y.*prob)));
    gOutput =  fp{end}.'*err;
    
    % If have more than one hidden layer, compute hidden layer weight gradients
    if length(nHidden) > 1
        
        % Last layer of hidden weights
        backprop = sech(ip{end}).^2.*(err*outputWeights');
        gHidden{end} = fp{end-1}'*backprop;
        
        % Other hidden layer weights
        for h = length(nHidden)-2:-1:1
            backprop =  sech(ip{h+1}).^2.*(backprop*hiddenWeights{h+1}');
            gHidden{h} =  fp{h}'*backprop;
        end
        
        % Input Weights
        backprop = sech(ip{1}).^2.*(backprop*hiddenWeights{1}');
        gInput   =  X'*backprop;
        
    elseif length(nHidden==1)
        
        % If have one hidden layer, compute only input weight gradients
        gInput =  X'*(sech(ip{end}).^2.*(err*outputWeights'));
    end
end

% Put gradient into vector
if nargout > 1
    g = zeros(size(w));
    if(~isempty(nHidden))
        g(1:nFeatures*nHidden(1)) = gInput(:);
        offset = nFeatures*nHidden(1);
        for h = 2:length(nHidden)
            g(offset+1:offset+nHidden(h-1)*nHidden(h)) = gHidden{h-1};
            offset = offset+nHidden(h-1)*nHidden(h);
        end
    else
        offset = 0;
    end
    g(offset+1:offset+nHidden(end)) = gOutput;
end
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
w = model.w;
nHidden = model.nHidden;

% Standardize optionally
if model.standardize
    Xhat = Xhat - repmat(model.mu, nTest, 1);
    Xhat = Xhat ./ repmat(model.sigma, nTest, 1);
end

% Extract weights from weight vector
if(~isempty(nHidden))
    inputWeights = reshape(w(1:nFeatures*nHidden(1)),nFeatures,nHidden(1));
    offset = nFeatures*nHidden(1);
    for h = 2:length(nHidden)
        hiddenWeights{h-1} = reshape(w(offset+1:offset+nHidden(h-1)*nHidden(h)),nHidden(h-1),nHidden(h));
        offset = offset+nHidden(h-1)*nHidden(h);
    end
else
    offset = 0;
end
outputWeights = w(offset+1:offset+nHidden(end));

% Compute Output
if(~isempty(nHidden))
    ip{1} = Xhat*inputWeights;
    % Hidden unit bias
    ip{1}(:,1) = -inf;
    fp{1} = tanh(ip{1});
    
    % Hidden unit bias
    fp{1}(:,1) = 1;
    for h = 2:length(nHidden)
        ip{h} = fp{h-1}*hiddenWeights{h-1};
        % Hidden unit bias
        ip{h}(:,1) = -inf;
        fp{h} = tanh(ip{h});
        % Hidden unit bias
        fp{h}(:,1) = 1;
    end
    prob  = sigmoid(fp{end}*outputWeights);
else
    prob  = sigmoid(Xhat*outputWeights{1});
end

% Final results
index = prob >= 0.5;
yhat = sign(index - 0.5);
end

function [phi] = sigmoid(x)
% Sigmoid function
phi = 1 ./ (1 + exp(-x));
end