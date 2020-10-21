function [model] = ml_unsupervised_sparseAutoencoder(X,options)
% ml_unsupervised_sparseAutoencoder(X,y,options)
%
% Description:
%	 Fit a feed-forward neural network with linear neurons using sigmoid
%	 activiation function, using reconstruction error (Frobenius norm) as
%	 loss function. Sparsity can be added to loss function through
%	 parameters, expressed as a penalty based on the KL-divergence of the
%	 hidden units with Bernoulli
%
% Options:
%   - nHidden: a vector of the size of each hidden layer,
%       including bias variable (default: [3, 3, 3])
%
%   - beta: sparsity parameter weight for each hidden layer. Set beta(i) to 
%           0 for no sparsity penalty at hidden layer i.
%   - rhos: average activitation of hidden units, modelled as Bernoulli
%           random variable. Set rho(i) to x in [0,1] to penalize layer i
%           for having a different activation than Bernouli(x) in the
%           KL-divergence sense.
%   - lambda: output weights penalty in L2 norm
%   - standardize: if set to 1, standardize columns of X to have
%                  mu 0, sigma 1 (default: 0)

[nTrain,nFeatures] = size(X);

% process inputs
[standardize, nHidden, lambda, sparsify, betas, rhos, gradcheck] = ...
    myProcessOptions(options, 'standardize', 0, 'nHidden', [196], ...
                     'lambda', 3e-3, 'sparsify', [1], 'betas', [3], ...
                     'rhos', [0.1], 'gradcheck', 0);

if ~isempty(sparsify) || ~isempty(betas) || ~isempty(rhos) 
    % validate sparsity parameters
    if length(nHidden)~=length(sparsify) || ...
       length(nHidden)~=length(betas) || ...
       length(nHidden)~= length(rhos)
        error(['Dimension mismatch for sparsity parameters: ', ...
            'sparsify must be logical vector same length as nHidden', ...
            '; rhos and betas vectors must be of same length as nHidden']);
    end
end

% Optimalization options
optimOptions.Display = 'off';
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

theta = initializeParameters(nFeatures, nHidden);

% Optimalization options
% Verify gradient implementation
if gradcheck
    % check gradient implementation 
    numgrad = computeNumericalGradient( ...
            @(T) sparseAutoencoderLoss(T, X,nHidden, lambda, sparsify, ...
                                       betas, rhos), theta);
    [~, grad] = sparseAutoencoderLoss(theta, X, nHidden, lambda, sparsify, ...
                                                betas, rhos);    
    diff = norm(numgrad-grad)/norm(numgrad+grad);
disp(diff);
end
optimOptions.Method = 'lbfgs'; 
optimOptions.maxIter = 500;
optimOptions.display = 'off';

% Optimize weights
w = minFunc(@sparseAutoencoderLoss,theta,optimOptions,X,nHidden,lambda, ...
            sparsify, betas, rhos);

% Model outputs
model.nHidden = nHidden;
model.name = 'Sparse Autoencoder';
model.w = w;
model.predict = @predict;
model.standardize = standardize;
end

function [f,g] = sparseAutoencoderLoss(w, X, nHidden, lambda, sparsify, ...
                                       betas, rhos)
% ----- INITIALIZE VARIABLES ----------------------------------------------

% Get number of instances, variables, and layers
[nTrain, nFeatures] = size(X);
L = 2+length(nHidden); % number of layers in network

[W, B] = extractParameters(w, nHidden, nFeatures);

% Initialize function and gradient
% For an autoencoder, there are L-1 gradient matrices for the L1\
for layer = 1:L-1
    gradW{layer} = zeros(size(W{layer}));
    gradB{layer} = zeros(size(W{layer},2));
end

% Initialize error deltas for backpropagation
D = cell(L,1);

% ----- FEED INPUTS FORWARDS ----------------------------------------------
% (Z{k})_{i,j} holds sum of inputs from layer k-1 to layer k from
% observation i to the j'th unit
% (A{k})_{i,j} holds the activations for observation i at the j'th unit
% at layer k

% Layer 1 is the input layer
A{1} = X;
for layer = 2:(L-1)
    Z{layer} = bsxfun(@plus, A{layer-1}*W{layer-1}, B{layer-1}');
    A{layer} = sigm(Z{layer});
end
Z{L} = bsxfun(@plus, A{L-1}*W{L-1}, B{L-1}');
A{L} = sigm(Z{L});

% ----- EVALUATE COST FUNCTION -------------------------------------------
relError = X - A{L};
f = (0.5)*(1/nTrain)*sum(relError(:).^2); % Frobenius matrix norm

% Add weight decay to cost function
for layer = 1:(L-1)
    f = f + 0.5*lambda*sum(W{layer}(:).^2);
end

% Encourage sparsity by KL divergence from Bern(rhos(i)) at hiddenLayer
for hiddenLayer = 1:length(nHidden);
    if sparsify(hiddenLayer)
        layerIdx = hiddenLayer + 1;
        Rhat{hiddenLayer} = sum(A{layerIdx}, 1)./nTrain; % Activations
        f = f + betas(hiddenLayer)*sum(KL(rhos(hiddenLayer), ...
                                          Rhat{hiddenLayer}));
    end
end

% ----- PROPAGATE ERROR DELTAS BACKWARDS ---------------------------------
D{L} = -relError.*ddxSigm(A{L});
for layer = L-1:-1:2 
    % adjust index for hidden layer data structures
    hiddenIdx = layer-1;
    if sparsify(hiddenIdx)
        %  activations are penalized at this layer
        ddrSparsity = betas(hiddenIdx)*ddrKL(rhos(hiddenIdx), ...
                                             Rhat{hiddenIdx});
        D{layer} = bsxfun(@plus, D{layer+1}*W{layer}', ...
                          ddrSparsity).*ddxSigm(A{layer});
    else
        % activations are not penalized at this layer
        D{layer} = D{layer+1}*W{layer}'.*ddxSigm(A{layer});
    end
end

% ----- EVALUATE GRADIENT ------------------------------------------------
% Evaluate gradient using backpropagated error deltas
gradW{L-1} = A{L-1}'*D{L}./nTrain + lambda*W{L-1};
gradB{L-1} = sum(D{L},1)./nTrain;
for layer = (L-2):-1:1
    gradW{layer} = A{layer}'*D{layer+1}./nTrain + lambda*W{layer};
    gradB{layer} = sum(D{layer+1},1)./nTrain;
end

g = unrollParameters(L, gradW, gradB);
end

function [theta] = initializeParameters(nFeatures, nHidden)
    % Initialize weights randomly within [-bound, bound]
    bound = sqrt(6/(nFeatures+nHidden(1) + 1));
    % Add input weights
    theta = 2*bound*rand(nFeatures*nHidden(1), 1) - bound;

    % Initialize weight matrices 
    for i = 2:length(nHidden);
        bound = sqrt(6/(nHidden(i-1)+nHidden(i) + 1));
        % Add hidden weights
        theta = [theta; 2*bound*rand(nHidden(i-1)*nHidden(i),1) - bound];
    end
    bound = sqrt(6/(nHidden(end)+nFeatures + 1));
    % Add output weights
    theta = [theta;  2*bound*rand(nHidden(end)*nFeatures,1) - bound];

    % Initialize bias terms
    for i = 1:length(nHidden)
        theta = [theta; zeros(nHidden(i),1)]; % Bias for hidden layers
    end
       theta = [theta; zeros(nFeatures,1)]; % Bias for output layer
end

function [U] = predict(model,Xhat)
% Prediction function
[~,nFeatures] = size(Xhat);
w = model.w;
nHidden = model.nHidden;

[W, B] = extractParameters(w, nHidden, nFeatures);
L = 2+length(nHidden);
% Compute Output
A{1} = Xhat;
for layer = 2:(L-1)
    Z{layer} = bsxfun(@plus, A{layer-1}*W{layer-1}, B{layer-1}');
    A{layer} = sigm(Z{layer});
end
Z{L} = bsxfun(@plus, A{L-1}*W{L-1}, B{L-1}');
U = sigm(Z{L});
end

%% Helper functions for sparse autoencoder
function [theta] = unrollParameters(numLayers, gradW, gradB)
% Unroll weight gradient
theta = gradW{1}(:);
for layer = 2:numLayers-1 
    theta = [theta; gradW{layer}(:)];
end
% Unroll bias gradient 
for layer = 1:(numLayers-1)
    theta = [theta; gradB{layer}(:)];
end
end

function [W, B] = extractParameters(theta, nHidden, nFeatures)
L = 2+length(nHidden);
% Extract weights from parameter vector
W{1} = reshape(theta(1:nFeatures*nHidden(1)), [nFeatures, nHidden(1)]); %nFeatures by n_{L-1}
offset = nFeatures*nHidden(1);
for l = 2:length(nHidden)
    W{l} = reshape(theta(offset+1:offset+nHidden(l-1)*nHidden(l)), ...
        [nHidden(l-1), nHidden(l)]); 
    offset = offset+nHidden(l-1)*nHidden(l);
end
W{L-1} = reshape(theta(offset+1:offset+nHidden(end)*nFeatures), ...
    [nHidden(end), nFeatures]); 
offset = offset + nHidden(end)*nFeatures;
% Extract bias from parameter vector
for l = 1:length(nHidden)
    B{l} = theta(offset+1:offset + nHidden(l));
    offset = offset + nHidden(l);
end
B{L-1} = theta(offset+1:offset +nFeatures); 
end

function [s] = sigm(X)
s =  1 ./ (1+exp(-X));
end

function [d] = ddxSigm(z)
d = z.*(1-z);
end

function [kl] = KL(r,Rhat) 
    kl = r.*log(r./Rhat)+(1-r).*log((1-r)./(1-Rhat));
end

function [d] = ddrKL(r, Rhat)
    d = -r./Rhat + (1-r)./(1-Rhat);
end