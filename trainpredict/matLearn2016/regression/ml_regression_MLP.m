function [model] = ml_regression_MLP(X,y,options)
% ml_regression_MLP(X,y,options)
%
% Description:
%   - Fits a multilayer perceptron model for regression
%
% Options:
%    - nHidden: a vector of the size of each hidden layer,
%       including bias variable (default: [3,3,3])
%    - activFunc: activation function used between each
%       layer, one of:
%       'tanh' - tanh activation function
%       'sig' - sigmoid activation function
%
% Authors:

[nTrain,nFeatures] = size(X);

% Default options
[nHidden, activFunc] = myProcessOptions(options,'nHidden', [3,3,3], 'activFunc', 'tanh');

% Initialize weights
bound = sqrt(6/(nFeatures+nHidden(1)));
if strcmp(activFunc, 'sig')
    bound = 4*bound;
end
w = 2*bound*rand(nFeatures*nHidden(1),1) - bound;
for i = 2:length(nHidden);
    bound = sqrt(6/(nHidden(i-1)+nHidden(i)));
    if strcmp(activFunc, 'sig')
        bound = 4*bound;
    end
    w = [w; 2*bound*rand(nHidden(i-1)*nHidden(i),1) - bound];
end
bound = sqrt(6/(nHidden(end)+1));
w = [w;  2*bound*rand(nHidden(end),1) - bound];

% Optimalization options
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Optimize weights
w = minFunc(@MLP,w,optimOptions,X,y,nHidden,activFunc);

% Model outputs
model.nHidden = nHidden;
model.w = w;
model.predict = @predict;
if strcmp(activFunc, 'tanh')
model.name = 'Multi-layer Perceptron with tanh Activation Function';
elseif strcmp(activFunc, 'sig')
model.name = 'Multi-layer Perceptron with Sigmoid/Logistic Activation Function';
end
model.activFunc = activFunc;
end

function [f,g] = MLP(w,X,y,nHidden,activFunc)
% Check for improper input
if(any(isnan(w))); f=inf; g=zeros(size(w)); return; end;

% Get number of instances, variables, and layers
[nTrain,nFeatures] = size(X);

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

% Initialize function and gradient
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
    if strcmp(activFunc, 'tanh')
        fp{1} = tanh(ip{1});
    elseif strcmp(activFunc, 'sig')
        fp{1} = sigmoid(ip{1});
    end
    [h, ~] = size(inputWeights);
    if (h ~= 1) %
        % Correct bias unit
        ip{1}(:, 1) = -inf;
        fp{1}(:, 1) = 1;
    end
    for h = 2:length(nHidden)
        ip{h} = fp{h-1}*hiddenWeights{h-1};
        if strcmp(activFunc, 'tanh')
            fp{h} = tanh(ip{h});
        elseif strcmp(activFunc, 'sig')
            fp{h} = sigmoid(ip{h});
        end
        [h, ~] = size(inputWeights);
        if (h ~= 1)
            %Correct for bias unit
            ip{h}(:,1) = -inf;
            fp{h}(:,1) = 1;
        end
    end
    yhat = fp{end}*outputWeights;
else
    yhat = X*outputWeights{1};
end

% Compute error matrix
relativeErr = yhat-y;
f = sum(relativeErr(:).^2); % should still be fine...

% Compute gradient
if nargout > 1
    err = 2*relativeErr;
    
    % Output weight gradient
    gOutput =  fp{end}'*err;
    
    % If have more than one hidden layer, compute hidden layer weight gradients
    if length(nHidden) > 1
        
        % Last layer of hidden weights
        if strcmp(activFunc, 'tanh')
            backprop = sech(ip{end}).^2.*(err*outputWeights');
        elseif strcmp(activFunc, 'sig')
            s = sigmoid(ip{end});
            backprop = s.*(1-s).*(err*outputWeights');
            
        end
        gHidden{end} = fp{end-1}'*backprop;
        
        % Other hidden layer weights
        for h = length(nHidden)-2:-1:1
            if strcmp(activFunc, 'tanh')
                backprop =  sech(ip{h+1}).^2.*(backprop*hiddenWeights{h+1}');
            elseif strcmp(activFunc, 'sig')
                s = sigmoid(ip{h+1});
                backprop =  s.*(1-s).*(backprop*hiddenWeights{h+1}');
            end
            gHidden{h} =  fp{h}'*backprop;
        end
        
        % Input Weights
        if strcmp(activFunc, 'tanh')
            backprop = sech(ip{1}).^2.*(backprop*hiddenWeights{1}');
        elseif strcmp(activFunc, 'sig')
            s = sigmoid(ip{1});
            backprop =  s.*(1-s).*(backprop*hiddenWeights{1}');
        end
        gInput   =  X'*backprop;
        
    elseif length(nHidden==1)
        
        % If have one hidden layer, compute only input weight gradients
        if strcmp(activFunc, 'tanh')
            gInput =  X'*(sech(ip{end}).^2.*(err*outputWeights'));
        elseif strcmp(activFunc, 'sig')
            s = sigmoid(ip{end});
            gInput =  X'*((s.*(1-s).*err)*outputWeights');
        end
        
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
    if strcmp(model.activFunc, 'tanh')
        fp{1} = tanh(ip{1});
    elseif strcmp(model.activFunc, 'sig')
        fp{1} = sigmoid(ip{1});
    end
    
    % Hidden unit bias
    fp{1}(:,1) = 1;
    for h = 2:length(nHidden)
        ip{h} = fp{h-1}*hiddenWeights{h-1};
        % Hidden unit bias
        ip{h}(:,1) = -inf;
        if strcmp(model.activFunc, 'tanh')
            fp{h} = tanh(ip{h});
        elseif strcmp(model.activFunc, 'sig')
            fp{h} = sigmoid(ip{h});
        end
        % Hidden unit bias
        fp{h}(:,1) = 1;
    end
    yhat  = fp{end}*outputWeights;
else
    yhat  = Xhat*outputWeights{1};
end
end

function [phi] = sigmoid(x)
% Sigmoid function
phi = 1 ./ (1 + exp(-x));
end
