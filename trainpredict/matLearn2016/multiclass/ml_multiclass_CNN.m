function [model] = ml_multiclass_CNN(X,y,options)
% ml_multiclass_CNN(X,y,options)
%
% Description:
% Find optimal set of parameters for simple convolutional neural network in
% multiclass classification regime.
%
% Options:
%   - imageDim:     Size of input image. Must be square. Default: 28.
%   - nClasses:     Number of classes in image. Default: 10.
%   - filterDim:    Size of (square) filters. Default: 9.
%   - nFilters:     Number of filters to use in convolutional layer.
%                   Default: 20.
%   - poolDim       Side length of square pooling region for mean pooling.
%                   Default: 20.
%   - convKind      Sets whether filter response should be 'valid', with
%                   no zero padding image edges, or 'same', with
%                   zero-padding. Defaut: 'valid'. (see conv2
%                   documentation)

[numClasses, filterDim, numFilters, poolDim, imageDim, convKind] = ...
            myProcessOptions(options, 'nClasses', -1, 'filterDim', -1, ...
                             'nFilters', -1, 'poolDim', -1, ...
                             'imageDim', -1, 'convKind', 'valid');

validateInput(numClasses, filterDim, numFilters, poolDim, imageDim);
% set convolution output dimension according to convolution type
if strcmp(convKind, 'same')
    outputDim = imageDim;
else
    outputDim = imageDim - filterDim + 1;
end
w = initializeParams(imageDim, filterDim, numFilters, poolDim, ...
                     numClasses, outputDim);              
X = reshape(X',imageDim,imageDim,[]); % reshape as tensor of images
nTrain= size(X,3);

% Stochastic Gradient Descent with momentum
maxIter = 100000;
stepSize = 1e-2;
funObj = @(w, i)cnnLoss(w, X(:,:,i), y(i), numClasses, filterDim, ...
                        numFilters, poolDim, false, convKind);
w_prev = w;
momentum = 0.95;
for iter = 1:maxIter
    if mod(iter - 1, round(maxIter / 10)) == 0
        fprintf('Training example: %d\n', iter);
    end
    i = ceil(rand*nTrain);
    [~, grad] = funObj(w,i);
    w_temp = w;
    w = w - stepSize*grad + momentum * (w - w_prev);    
    w_prev = w_temp;
end                   

model.w = w;
model.convKind = convKind;
model.imageDim = imageDim;
model.nClasses = numClasses;
model.filterDim = filterDim;
model.nFilters = numFilters;
model.poolDim = poolDim;
model.predict = @predict;
model.name = 'CNN for Multiclass Image Classification';
end

function [preds] = predict(model, Xtest)
Xtest = reshape(Xtest',model.imageDim,model.imageDim,[]);
[~,~,preds] = cnnLoss(model.w, Xtest, [], model.nClasses, model.filterDim, ...
                model.nFilters, model.poolDim, true, model.convKind);
end

function [f, grad, preds] = cnnLoss(w, images, labels, nClasses, ...
                                    filterDim, nFilters, poolDim, ...
                                    pred, convKind)
inputImageDim = size(images,1); % height / width of one square image
nImages = size(images,3);

%% ----- FORWARD PASS -----------------------------------------------------
% Convolutional layer -----------------------------------------------------
if strcmp(convKind, 'same') % treat image as zero-padded for convolution?
    convDim = inputImageDim; % dimension of convolved output
else
    convDim = inputImageDim - filterDim + 1;
end

outputDim = convDim / poolDim; % dimension of subsampled output
[Wc, Wd, bc, bd] = reshapeParams(w, outputDim, filterDim, ...
                                         nFilters,nClasses);

% convDim x convDim x numFilters x numImages tensor of filter activations
filterActivations = convolve(filterDim, nFilters, images, Wc, bc, convKind);

% outputDim x outputDim x numFilters x numImages tensor of
% subsampled activations
pooledActivations = pool(poolDim, filterActivations);

% Pooling layer -----------------------------------------------------------
% reshape activations into 2-d matrix, hiddenSize x numImages,
% for Softmax layer
pooledActivations = reshape(pooledActivations, [], nImages);

% Softmax Layer -----------------------------------------------------------
% denominator: row vector of 1 x nImages containing the sum for each class
unnormedLogProb = bsxfun(@plus, Wd*pooledActivations, bd);

% use log-sum-exp computational trick to prevent overflow/underflow
normedLogProbs = bsxfun(@minus, unnormedLogProb, ...
                        logSumExp(unnormedLogProb));

% reuse forward pass code for model's predict function
if pred
    [~,preds] = max(exp(normedLogProbs),[],1);
    preds = preds';
    grad = 0;
    f = 0;
    return;
end

% Evaluate cost function --------------------------------------------------
onesIdx = sub2ind(size(unnormedLogProb), labels', 1:nImages);
trainingDist = zeros(size(unnormedLogProb));
trainingDist(onesIdx) = 1;

% average cross-entropy loss over training examples
f = -(1/nImages)*sum(sum(normedLogProbs .* trainingDist, 1));

% derivative of loss w.r.t. input to softmax layer
softmaxErr = (1/nImages)*(exp(normedLogProbs) - trainingDist);

%% ----- PROPAGATE ERRORS BACKWARDS ---------------------------------------
% gradient of softmax loss w.r.t. output of pooling layer
poolingOutputErr = Wd' * softmaxErr;

% reshape as 4D array for upsample operation
poolingOutputErr = reshape(poolingOutputErr, [], outputDim, nFilters, ...
                           nImages);  
                       
% gradient of loss w.r.t input of pooling layer                         
poolingInputErr = upsampleMeanPoolErr(poolingOutputErr, poolDim , ...
                                      convDim, nImages, nFilters); 
                                     
% gradient of loss w.r.t output of convolution, prior to sigmoid activation
convOutErr = poolingInputErr .* ddx_sigm(filterActivations);

%% ----- FIND GRADIENT W.R.T. PARAMETERS ----------------------------------
% grad of loss w.r.t. weights and biases of softmax input
Wd_grad = softmaxErr * pooledActivations.';
bd_grad = sum(softmaxErr, 2);

% grad of loss w.r.t. filters and biases of convolution input
Wc_grad = backpropConvErr(convOutErr, images, nImages, filterDim, ...
                          nFilters);
bc_grad = squeeze(sum(sum(sum(convOutErr, 4), 2), 1));

grad = [Wc_grad(:); Wd_grad(:); bc_grad(:); bd_grad(:)];
end

function[Wc, Wd, bc, bd] = reshapeParams(w, convOutDim, filterDim, ...
                                         numFilters,numClasses)
hiddenSize = convOutDim^2*numFilters;

% Read filters from flat list and reshape
offset = filterDim^2*numFilters;
Wc = reshape(w(1:offset),filterDim,filterDim,numFilters);
Wd = reshape(w(offset + 1:offset + hiddenSize*numClasses),numClasses,hiddenSize);

% Read bias vectors from flat list
offset = offset + hiddenSize*numClasses;
bc = w(offset + 1:offset + numFilters);
offset = offset + numFilters;
bd = w(offset + 1:offset + numClasses);
end

function[theta] = initializeParams(imageDim, filterDim, numFilters, ...
    poolDim, numClasses, outputDim)
% TODO check: filterDim is less than imageDim
% initializeParameters
Wc = randn(filterDim, filterDim, numFilters)./sqrt(filterDim.^2*numFilters);
% TODO check: outDim is a multiple of poolDim
outputDim = outputDim/poolDim;
hiddenSize = outputDim^2*numFilters; % number of parameters in the hidden layer
bound = sqrt(6/(numClasses+hiddenSize + 1));
Wd = 2*bound*rand(numClasses,hiddenSize) - bound;
bc = zeros(numFilters, 1);
bd = zeros(numClasses, 1);
theta = [Wc(:); Wd(:); bc(:); bd(:)];
end

function [s] = sigm(X)
s = 1 ./ (1 + exp(-X));
end

function [d] = ddx_sigm(X)
d = X .* (1 - X);
end

function [convInputErr] = backpropConvErr(convOutErr, images, numImages, ...
                                          filterDim, numFilters)
convInputErr = zeros(filterDim, filterDim, numFilters);
for f = 1:numFilters
    gradAcc = zeros(filterDim);
    for i = 1:numImages
        % convolve each image with gradient from output layer to sum up 
        % contributions to gradient from each image 
        flippedErrFilter = rot90(convOutErr(:,:,f,i), 2);
        gradAcc = gradAcc+conv2(images(:,:,i),flippedErrFilter,'valid');
    end
    convInputErr(:,:,f) = gradAcc;
end
end

function [err] = upsampleMeanPoolErr(outputErrors, poolDim, convDim, ...
                                        numImages, numFilters)                                  
err = zeros(convDim, convDim, numFilters, numImages);
unpoolFilter = ones(poolDim,poolDim) / (poolDim^2);
for i = 1:numImages
    for f = 1:numFilters
        % upsample errors through mean pooling layer
        err(:,:,f,i) = kron(outputErrors(:,:,f,i), unpoolFilter);
    end
end                                                          
end

function [lse] = logSumExp(X)
% computes log sum exponential over the rows of a matrix X
Xmax = max(X);
lse = log(sum(exp(bsxfun(@minus, X, Xmax)),1)) + Xmax;
end

function [pooledFeatures] = pool(poolDim, convolvedFeatures)
% REQUIRES: filter is square and <= to the feature size
% find parameters of pooling operation
numImages = size(convolvedFeatures, 4);
numFilters = size(convolvedFeatures, 3);
convolvedDim = size(convolvedFeatures, 1);
filterDim = convolvedDim / poolDim;
% we only need the pooled features for non-overlapping regions, at the 
% indices in the "targets" vector:
targets = 1:poolDim:convolvedDim-1;

pooledFeatures = zeros(filterDim, filterDim, numFilters, numImages);
meanFilter = ones(poolDim, poolDim) / poolDim^2;
for i = 1:numImages
    for f = 1:numFilters                
      feature = convolvedFeatures(:, :, f, i);
      convolvedFeature = conv2(feature, meanFilter, 'valid');
      pooledFeatures(:, :, f, i) = convolvedFeature(targets, targets);             
    end
end
end

function [convolvedFeatures] = convolve(filterDim, numFilters, images, ...
                                        W, b, convKind)
numImages = size(images, 3);
imageDim = size(images, 1);
if strcmp(convKind, 'same')
    convDim = imageDim;
else
    convDim = imageDim - filterDim + 1;
end
convolvedFeatures = zeros(convDim, convDim, numFilters, numImages);

for imageNum = 1:numImages
  for filterNum = 1:numFilters
    im = squeeze(images(:, :, imageNum));
    filter = squeeze(W(:,:,filterNum));
        convolvedImage = zeros(convDim, convDim);%zeros(size(im));
    
    % a discrete convolution by definition calculates the dot product as 
    % if the image were flipped horizontally and vertically, so here we
    % correct that here in order to output the correct filter response
    filter = rot90(filter,2);
    convolvedImage = bsxfun(@plus, convolvedImage, conv2(im, filter, ...
                                                         convKind));          
    convolvedImage = sigm(bsxfun(@plus, convolvedImage, b(filterNum)));
    convolvedFeatures(:, :, filterNum, imageNum) = convolvedImage;
  end
end

end
function [] = errorMsg(str)
error('Error: must specify %s', str);
end

function [] = validateInput(numClasses, filterDim, numFilters, poolDim, ...
                            imageDim)
errorFlag = false;
if imageDim == -1
    errorMsg('dimension of (square) image');
    errorFlag = true;
end;
if numClasses == -1
    errorMsg('number of classes');
    errorFlag = true;
end;
if filterDim == -1
    errorMsg('number of filters');
    errorFlag = true;
end
if numFilters== -1
    errorMsg('filter dimension');
    errorFlag = true;
end
if poolDim == -1
    errorMsg('pooling dimension');
    errorFlag = true;
end
if errorFlag
    return;
end
end
