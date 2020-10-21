function [model] = ml_multilabel_MLP(X,y,options)
% ml_multilabel_MLP(X,y,options)
%
% Description:
%	 - Multilabel classification using Multi-Layer Perceptrons (neural
%      networks) 
%
% Options:
%   - nLabels: number of labels in dataset (default: None)
%   - verbose: 1 for detailed optimizer output, 0 no output (default: 1)
%   - nHidden: number of entries in vector gives number of hidden layers,
%              quantity of each entry gives number of hidden units in layer
%              (default: 3)
%   - lambdaI: L2-regularization for input layer weights (default: 0)
%   - lambdaH: L2-regularization for hidden layer weights (default: 0)
%   - lambdaO: L2-regularization for output layer (default: 0)
%
% Authors:
% 	 - Mark Schmidt (2014); adapted for matLearn by Geoffrey Roeder (2016)

if nargin < 3
    options = [];
end

[nLabels,verbose,nHidden,lambdaI,lambdaH,lambdaO] = ...
    myProcessOptions(options,'nLabels',[],'verbose',0,'nHidden',3,...
                     'lambdaI',0,'lambdaH',0,'lambdaO',0);

if ~verbose
    options.Display = 'none';
end

if isempty(nLabels)
    assert(1==0,'nLabels must be specified for multi-label classifiers\n');
end

y = linearInd2Binary(y,nLabels);

[nInstances,nVars] = size(X);

loss = @(w)ml_MLP_multiple_regression_loss(w,X,y,nHidden,nLabels);

nParams = nVars*nHidden(1);
for h = 2:length(nHidden);
    nParams = nParams+nHidden(h-1)*nHidden(h);
end
nParams = nParams+nHidden(end)*nLabels;

w_init = randn(nParams,1);
%options.numDiff = 1;
%options.DerivativeCheck = 'on';
if any([lambdaI;lambdaH;lambdaO] > 0)
    % Set up regularizer
    lambda(1:nVars*nHidden(1),1) = lambdaI;
    offset = nVars*nHidden(1);
    for h = 2:length(nHidden)
        lambda(offset+1:offset+nHidden(h-1)*nHidden(h),1) = lambdaH;
        offset = offset+nHidden(h-1)*nHidden(h);
    end
    lambda(offset+1:offset+nHidden(end)*nLabels,1) = lambdaO;
    penalizedLoss = @(w)ml_penalized_L2(w,loss,lambda);
    w = minFunc(penalizedLoss,w_init,options);
else
    w = minFunc(loss,w_init,options);
end

model.name = 'Multi-Label MLP';
model.nClasses = 2^nLabels;
model.nLabels = nLabels;
model.weights = w;
model.nHidden = nHidden;
model.predict = @predict;
end

function [y] = predict(model,X)
nHidden = model.nHidden;
w = model.weights;
nLabels = model.nLabels;

[nInstances,nVars] = size(X);

% Form Weights
inputWeights = reshape(w(1:nVars*nHidden(1)),nVars,nHidden(1));
offset = nVars*nHidden(1);
for h = 2:length(nHidden)
  hiddenWeights{h-1} = reshape(w(offset+1:offset+nHidden(h-1)*nHidden(h)),nHidden(h-1),nHidden(h));
  offset = offset+nHidden(h-1)*nHidden(h);
end
outputWeights = w(offset+1:offset+nHidden(end)*nLabels);
outputWeights = reshape(outputWeights,nHidden(end),nLabels);

% Compute Output
for i = 1:nInstances
    ip{1} = X(i,:)*inputWeights;
    fp{1} = tanh(ip{1});
    for h = 2:length(nHidden)
        ip{h} = fp{h-1}*hiddenWeights{h-1};
        fp{h} = tanh(ip{h});
    end
    y(i,:) = sign(fp{end}*outputWeights);
end
y = binary2LinearInd(y);
end