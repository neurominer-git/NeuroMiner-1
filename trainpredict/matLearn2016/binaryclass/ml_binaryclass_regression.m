function [model] = ml_binaryclass_regression(X,y,options)
% ml_binaryclass_regression(X,y,options)
%
% Description:
%	 - Classification by using a regression model
%
% Options:
%    - subModel: submodel on which the classifier needs to be trained
%       (default: mean)
%	 - getScores: whether to get score for each test containing the output
%       value obtained by using the submodel regressor (default: 0)
%	 - subOptions: options concerning the submodel (default: none)
%	 - crossValidate: specify number of cross validation bins to use
%       (default: 1)
%	 - bagging: specify whether to perform bagging on the training set
%       (default: 0)
%
% Authors:
% 	 - Anurag Ranjan (2014)

[nTrain,nFeatures] = size(X);

% Define options if not specified by input
[getScores,subModel,subOptions,crossValidate,bagging] = ... 
    myProcessOptions(options, 'getScores',0, 'subModel', ...
                     @ml_regression_mean,'subOptions',[], ...
                     'crossValidate',1,'bagging',0);

% Bagging (randomly pick and reorder training sets)
if bagging
    for i = 1:nTrain
        k = randi(nTrain);
        Xbag(i,:) = X(k,:);
        ybag(i,:) = y(k,:);
    end
    X = Xbag;
    y = ybag;
end

if (crossValidate>1)
    % Choose best regression model if crossValidate > 1
    for i = 1:crossValidate
        
        % Determine the indices for cross validation set
        validInd = (i-1)*(floor(nTrain/crossValidate))+1 ...
                   : (i)*(floor(nTrain/crossValidate));
        
        % Prepare test data (indexed data)
        XtestValid = X(validInd,:);
        ytestValid = y(validInd,:);
        
        % Prepare training data (all data - indexed)
        XtrainValid = X;
        XtrainValid(validInd,:) = [];
        ytrainValid = y;
        ytrainValid(validInd,:) = [];
        % Prepare submodel
        subModelValid{i} = subModel(XtrainValid,ytrainValid,subOptions);
        
        % Compute scores
        yhatValid = sign(subModelValid{i}.predict(subModelValid{i}, ...
                         XtestValid));
        testError(i) = sum(yhatValid~=ytestValid)/length(ytestValid);
    end
    
    % Find minimum validation error and best model
    [minVal,modelInd] = min(testError);
    regModel = subModelValid{modelInd};
else
    % Otherwise use use regression model once
    regModel = subModel(X,y,subOptions);
end

% Model outputs
model.regModel = regModel;
model.name = ['Binary Classification by ', regModel.name];
model.getScores = getScores;
model.predict = @predict;
model.crossValidate = crossValidate;
model.bagging = bagging;
model.submodelOptions = subOptions;
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);
yhat = zeros(nTest,1);

% Use regression predict function
sub_yhat = model.regModel.predict(model.regModel, Xhat);

index = sub_yhat >= 0;
yhat = sign(index - 0.5);
end