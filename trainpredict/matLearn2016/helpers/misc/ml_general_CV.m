function [model] = ml_general_CV(X, y, options)
% Description:
% 	 - This computes the "best" hyper-parameter(s) using cross-validation for
%       classification and regression problems
%       (e.g. nHidden with multi-layer perceptron)
%
% Options:
%	 - subModel: the learning algorithm (default: mean)
%	 - subOptions: options/parameters you want to pass through to the
%       learning algorithm (default: 0)
%	 - nFolds: number of folds (default: 5)
%	 - paramNames: name(s) of parameter(s) to optimize over
%       * single hyper-parameter - must be a single string
%       * two hyper-parameter - must be a cell array containing
%       the name of each paramater as a string e.g. {'lambda' , 'sigma'}
%       (default: none)
%    - paramValues: list of values to grid search against
%       * for single parameter CV, paramValues must be a vector
%       e.g. [1, 10, 100, 1000]
%       * for two paramater CV, paramValues must be a cell array of
%       size [1 2], where each column contains values corresponding to
%       each parameter in paramNames
%       e.g. lambdaValues = 2.^(2:-1:-12)';
%       sigmaValues = 2.^(3:-1:-7)';
%       paramValues = {lambdaValues , sigmaValues}
%       (default: none)
%    - loss: name of the loss function used to compute the
%       cross-validation score, one of:
%       'mn' - absolute error
%       'sq' - square error
%       'zo' - zero-one loss
%	'mc' - misclassification rate
%       (must be specified)
%	 - earlyStop: stop grid search when the error starts increasing after
%       it has decreased at least once (default: 0)
%	 - shuffle: whether to shuffle the data (default: 0)
%	 - nParams: whether the cross-validation is performed to
%       optimize over a single parameter (1) or two parameters (2)
%       (default: 1)
%	 - leaveOneOut: each sample is used once as a validation set
%       (singleton) while the remaining samples form the
%       training set; equivalent to setting nFolds to the number of
%       samples in the data (default: 0)
%    * loss, paramNames, and paramValues are mandatory options
%
% Output:
%	 - bestModel: the model with the parameter that achieved the best
%       score
%	 - bestParamValue: the best parameter value found
%    - bestError: the average validation error achieved with the best
%       parameter value
%	 - validationErrorLog: struct containing:
%       * paramValues: the parameter values list as used by CV (sorted)
%       * errorValues: the average validation error corresponding to the
%       parameter values
%
% Authors:
%	 - Issam Laradji (2014): issam.laradji@gmail.com
%	 - Matthew Dirks (2014): http://www.cs.ubc.ca/~mcdirks/
%	 - Nazanin Hamzei (2014): nazaninh@ece.ubc.ca

X_copy = X;
y_copy = y;
model.errorLog = [];
model.predict = @predict;

%     Default values
[subModel, subOptions, nFolds, paramNames, paramValues, ...
    loss, shuffle, earlyStop, nParams, leaveOneOut] ...
    = myProcessOptions(options,  ...
    'subModel', @ml_regression_mean,               ...
    'subOptions', [],         ...
    'nFolds', 2,                ...
    'paramNames', NaN,           ...
    'paramValues', NaN,         ...
    'loss', NaN,     ...
    'shuffle', 0,           ...
    'earlyStop', 0,         ...
    'nParams',1,         ...
    'leaveOneOut', 0        );
% Default return values
model.subModel = NaN;
if nParams == 1;
    model.name = sprintf('CV on: \r %s', paramNames);
else
    model.name = sprintf('CV on: %s & %s',paramNames{1},paramNames{2});
end
model.bestParams = NaN;
model.minError = NaN;
% Verify mandatory arguments
if ~isa(subModel,'function_handle')
    fprintf('ERROR: model must be specified, and must be a function handle.\n');
    return;
end

if isnan(loss)
    fprintf('ERROR: loss must be specified. If unsure, use ''sq'' for regression or ''mc'' for classification.\n')
end
if nParams==1 & isnan(paramNames)
    fprintf('ERROR: paramNames is a mandatory option which was not specified.\n');
    return;
end
if nParams==1 & ~ischar(paramNames)
    fprintf('ERROR: paramNames must be a single string.\n');
    return;
end
if nParams==1 & isscalar(paramValues)
    if isnan(paramValues)
        fprintf('ERROR: paramValues is a mandatory option which was not specified.\n');
        return;
    end
    fprintf('ERROR: paramValues must be a vector of more than 1 value.\n');
    return;
end
if nParams==2 & ~iscell(paramNames)
    fprintf('ERROR: paramNames must be a cell array.\n')
    return;
end
%     if max(isnan(cell2mat(paramNames)))
if nParams==2 & max(cell2mat(cellfun(@isnan,paramNames,'UniformOutput',false)))
    fprintf('ERROR: paramNames is a mandatory option which was not specified.\n');
    return;
end

if nParams==2
    for checkParamNames = 1:size(paramNames,2)
        if ~ischar(paramNames{checkParamNames})
            fprintf('ERROR: each entry in paramNames must be a string.\n')
            return;
        end
    end
end

if nParams==2
    for checkParamValues = 1:size(paramValues,2)
        if isscalar(paramValues{checkParamValues})
            if isnan(paramValues{checkParamValues})
                fprintf('ERROR: paramValues is a mandatory option which was not specified.\n');
                return;
            end
            fprintf('ERROR: each element(column) of paramValues must be a vector of more than 1 value.\n');
            return;
        end
    end
end

% Verify valid nFolds option
if nFolds > size(X,1)
    fprintf('WARNING: nFolds (%d) must not exceed number of samples in the data (%d), using %d for nFolds instead.\n', nFolds, size(X,1), size(X,1));
    nFolds = size(X,1);
end
if nFolds < 2
    fprintf('WARNING: nFolds (%d) must be greater than 1, using 2 for nFolds instead.\n', nFolds);
    nFolds = 2;
end

% Leave-one-out is a special case of k-fold CV
% where nFolds = number of samples
if leaveOneOut
    nFolds = size(X,1);
end

% Set loss function
if strcmp(loss, 'sq')
    lossFunction = @square_error;
elseif strcmp(loss, 'mc')
    lossFunction = @misclass_rate;
elseif strcmp(loss, 'zo')
    lossFunction = @zero_one_loss;
elseif strcmp(loss, 'mn')
    lossFunction = @abs_error;
else % Default if user provides invalid option
    fprintf('WARNING: invalid loss function name.\n');
    lossFunction = @square_error;
end

% Shuffle dataset
if shuffle
    randIndices = randperm(length(X));
    X = X(randIndices, :);
    y = y(randIndices, :);
end

% Make sure the param values are sorted
if nParams==1
    paramValues = sort(paramValues);
else
    paramValues{1} = sort(paramValues{1});
    paramValues{2} = sort(paramValues{2});
end

model.minError = inf;
model.bestParams = NaN;
isDecreasedPrev = false;
model.errorLog.errorValues = [];
if nParams == 1;
    for paramIndex = 1:length(paramValues)
        subOptions = SetProp(subOptions,options.paramNames,paramValues(paramIndex));
        validationErrors = zeros(nFolds,1);
        
        % Compute the accumulated score over the folds
        for fold = 1:nFolds
            % Split dataset into training and validation sets
            [Xtrain, ytrain, Xval, yval] = foldData(X, y, nFolds, fold);
            
            % Train model
            trainedModel = subModel(Xtrain, ytrain, subOptions);
            
            % Predict y
            yhat = trainedModel.predict(trainedModel, Xval);
            
            % Record error for this fold
            validationErrors(fold) = lossFunction(yhat, yval);
        end
        
        avgValidationError = mean(validationErrors);
        model.errorLog.errorValues(paramIndex) = avgValidationError;
        
        % Naive early stop method to prune parameter search values
        if earlyStop && paramIndex > 1
            % Check if error is decreasing
            if prevValidationError > avgValidationError
                isDecreasedPrev = true;
            end
            % Check if error is increasing,
            % having decreased at least once previously
            if prevValidationError < avgValidationError && isDecreasedPrev
                break;
            end
        end
        prevValidationError = avgValidationError;
        
        % Check if new validation error beats previous best
        if avgValidationError < model.minError
            model.minError = avgValidationError;
            model.subModel = trainedModel;
            
            if iscell(paramValues)
                model.bestParams = paramValues{paramIndex};
            else
                model.bestParams = paramValues(paramIndex);
            end
        end
    end
else
    
    for param1index = 1:length(paramValues{1})
        for param2index = 1:length(paramValues{2})
            subOptions = SetProp(subOptions,options.paramNames{1},paramValues{1}(param1index));
            subOptions = SetProp(subOptions,options.paramNames{2},paramValues{2}(param2index));
            validationErrors = zeros(nFolds,1);
            
            % Compute the accumulated score over the folds
            for fold = 1:nFolds
                % Split dataset into training and validation sets
                [Xtrain, ytrain, Xval, yval] = foldData(X, y, nFolds, fold);
                
                % Train model
                trainedModel = subModel(Xtrain, ytrain, subOptions);
                
                % Predict y
                yhat = trainedModel.predict(trainedModel, Xval);
                
                % Record error for this fold
                validationErrors(fold) = lossFunction(yhat, yval);
            end
            
            avgValidationError = mean(validationErrors);
            model.errorLog.errorValues(param1index,param2index) = avgValidationError;
            
            % Naive early stop method to prune parameter search values
            if earlyStop && param1index > 1 && param2index > 1
                % Check if error is decreasing
                if prevValidationError > avgValidationError
                    isDecreasedPrev = true;
                end
                % Check if error is increasing,
                % having decreased at least once previously
                if prevValidationError < avgValidationError && isDecreasedPrev
                    break;
                end
            end
            prevValidationError = avgValidationError;
            
            % Check if new validation error beats previous best
            if avgValidationError < model.minError
                model.minError = avgValidationError;
                model.subModel = trainedModel;
                
                model.bestParams = [paramValues{1}(param1index),paramValues{2}(param2index)];
            end
        end
    end
end

if nParams==1
    subOptions = SetProp(subOptions,options.paramNames,model.bestParams);
else
    subOptions = SetProp(subOptions,options.paramNames{1},model.bestParams(1));
    subOptions = SetProp(subOptions,options.paramNames{2},model.bestParams(2));
    
    
end
model.subModel = subModel(X_copy, y_copy, subOptions);

% Return the set of parameters that were actually tried, and the order they were used.
if nParams==1
    model.errorLog.paramValues = paramValues(1:length(model.errorLog.errorValues));
else
    model.errorLog.paramValues.paramNames{1} = paramValues{1}(1:length(paramValues{1}));
    model.errorLog.paramValues.paramNames{2} = paramValues{1}(1:length(paramValues{1}));
end

end

function [Xtrain, ytrain, Xval, yval] = foldData(X, y, nFolds, fold)
% nFolds: total number of folds; fold: which fold to perform

n = size(y,1);

% Calculate number of samples per fold
% (except the last fold will take all remaining rows)
nPerFold = floor(n/nFolds);

% Create validation set
valStart = (fold-1)*nPerFold + 1;
if (fold == nFolds)
    valEnd = n;
else
    valEnd = fold*nPerFold;
end

yval = y(valStart:valEnd, :);
Xval = X(valStart:valEnd, :);

% Training set may consist of 2 parts: the part before the validation
% set, and the part after
if (fold == 1)
    trainStart = valEnd + 1;
    trainEnd = n;
    ytrain = y(trainStart:trainEnd, :);
    Xtrain = X(trainStart:trainEnd, :);
elseif (fold == nFolds)
    trainStart = 1;
    trainEnd = valStart - 1;
    ytrain = y(trainStart:trainEnd, :);
    Xtrain = X(trainStart:trainEnd, :);
else
    % Expect 2 parts, because fold must be in the middle
    trainStartA = 1;
    trainEndA = valStart - 1;
    trainStartB = valEnd + 1;
    trainEndB = n;
    
    ytrain = [
        y(trainStartA:trainEndA, :);
        y(trainStartB:trainEndB, :)
        ];
    Xtrain = [
        X(trainStartA:trainEndA, :);
        X(trainStartB:trainEndB, :)
        ];
end
end

function [error] = square_error(yhat, y)
% Also known as: average L2^2, or MSE
error = mean(mean((yhat - y).^2));
end

function [error] = abs_error(yhat, y)
% Also known as: average L1 norm
error = mean(mean(abs(yhat - y)));
end

function [error] = zero_one_loss(yhat, y)
error = mean(mean((yhat ~= y)));
end

function [error] = misclass_rate(yhat, y)
error = mean(mean((yhat ~= y)));
end


function [yhat] = predict(model, Xhat)
% Predict Function
yhat = model.subModel.predict(model.subModel, Xhat);
end

function [obj] = SetProp(obj, prop, std_val)
% Function that sets nested fields
S=struct('type','.','subs',regexp(prop,'\.','split'));
obj = subsasgn(obj, S, std_val);
end
