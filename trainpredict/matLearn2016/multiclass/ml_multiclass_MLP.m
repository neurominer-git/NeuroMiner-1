function [model] = ml_multiclass_MLP(X,y,options)
% ml_multiclass_MLP(X,y,options)
%
% Description:
%	 - Classification using a multilayer perceptron with softmax
%       loss
%
% Options:
%    - nHidden: a vector of the size of each hidden layer,
%       including bias variable (default: [3,3,3])
%    - standardize: if set to 1, standardize columns of X to have
%       mu 0, sigma 1 (default: 0)

[nTrain,nFeatures] = size(X);

% process inputs
[standardize,nHidden] =...
    myProcessOptions(options,...
    'standardize', 0, 'nHidden', 3*ones(1,3));

classes = unique(y);
nClasses = length(classes);

% Optimalization options
optimOptions.Display = 0;
optimOptions.useMex = 0;

if standardize
    mu = mean(X);
    sigma = std(X);
    X = X - repmat(mu, nTrain, 1);
    X = X ./ repmat(sigma, nTrain, 1);
    
    % Make sure we can apply same transformation at testing time
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
bound = sqrt(6/(nHidden(end)+nClasses));
w = [w;  2*bound*rand(nHidden(end)*nClasses,1) - bound];

% Optimalization options
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Optimize weights
w = minFunc(@MLP,w,optimOptions,X,y,nHidden, classes);

% Model outputs
model.nHidden = nHidden;
model.name = 'Multi-layer Perceptron Classification';
model.w = w;
model.predict = @predict;
model.standardize = standardize;
model.nClasses = nClasses;
model.classes = classes;
end

function [f,g] = MLP(w, X, y, nHidden, classes)
% Loss function for MLP

% Check for improper input
if(any(isnan(w))); f=inf; g=zeros(size(w)); return; end;

% Get number of instances, variables, and layers
[nTrain,nFeatures] = size(X);
nClasses = length(classes);

% Extract weights from weight vector
if (~isempty(nHidden)) % input 2 x 3: maps R^3 to R^2
    inputWeights = reshape(w(1:nFeatures*nHidden(1)),nFeatures,nHidden(1));
    offset = nFeatures*nHidden(1);
    for h = 2:length(nHidden) % 3 x 5 for one hidden layer: maps R^5 to R^3
        hiddenWeights{h-1} = reshape(w(offset+1:offset+nHidden(h-1)*nHidden(h)),nHidden(h-1),nHidden(h));
        offset = offset+nHidden(h-1)*nHidden(h);
    end
else
    offset = 0;
end
% maps R^5 to R^5
outputWeights = reshape(w(offset+1:offset+nHidden(end)*nClasses),nHidden(end),nClasses);

% Initialize function and gradiant
if nargout > 1
    if(~isempty(nHidden))
        gInput = zeros(size(inputWeights)); % initialize gradient in to zeros
    else
        gInput = [];
    end
    gOutput = zeros(size(outputWeights)); % intialize output layer gradient to zero (deltas)
    for h = 1:length(nHidden)-1
        gHidden{h} = zeros(size(hiddenWeights{h})); % initialize two hidden layers to zero
    end
end

% Compute Output
if(~isempty(nHidden))
    ip{1} = X*inputWeights; % ip holds z_j&(l) for layer l z_j^(1) = X * W_ij^(1) from 2 to 3
    fp{1} = tanh(ip{1}); % fp holds activiations a_1 
    
    % Correct bias unit
    ip{1}(:,1) = -inf; % what does the correction do here?
    fp{1}(:,1) = 1; % set the last value to 1
    for h = 2:length(nHidden)
        ip{h} = fp{h-1}*hiddenWeights{h-1}; % for each of the (2) hidden layers find z_i^(h) using a_i^(h-1)
        fp{h} = tanh(ip{h}); % find activation 
        
        %Correct for bias unit
        ip{h}(:,1) = -inf; % set the previous bias to -inf
        fp{h}(:,1) = 1; % lock the bias at this hidden layer to 1
    end
    prob = fp{end}*outputWeights; % find the final output 
else
    prob = X*outputWeights{1}; % there were no hidden layers, so find z_i
end

yBin = bsxfun(@eq, y, classes'); % generate a 0-1 logical vector for y
% prob - log exp outputweights^T * activation^(n_l)
tmp = prob - repmat(logsumexp(prob), 1, nClasses); % tmp is the second terms of the losf
% function J(theta). tmp hols for each example i and class k the log(.)
% term in the loss function
f = -sum(sum(tmp .* yBin, 2)); % find output using softmax regression

% Compute gradient
if nargout > 1
    
    % Output weight gradient
    mu = exp(tmp); % mu is the exp / sum (exp) term
    err =  -(yBin - mu);
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
    g(offset+1:offset+nHidden(end)*nClasses) = reshape(gOutput,nHidden(end) * nClasses, 1);
end
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
w = model.w;
nHidden = model.nHidden;

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
outputWeights =  reshape(w(offset+1:offset+nHidden(end)*model.nClasses),nHidden(end),model.nClasses);

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
    num  = fp{end}*outputWeights;
else
    num  = Xhat*outputWeights{1};
end

% Final results
[~, indx] = max(num, [], 2);
yhat = model.classes(indx);
end

function [lse] = logsumexp(b)
% Compute log(sum(exp)) across columns without overflowing
B = max(b,[],2);
lse = log(sum(exp(b-repmat(B,[1 size(b,2)])),2))+B;
end
