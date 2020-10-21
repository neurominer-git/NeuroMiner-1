function [model] = ml_multiclass_1v1(X,y,options)
% ml_multiclass_1v1(X,y,options)
%
% Description:
%    - Train C(C-1)/2 binary classifiers for a C-way multiclass problem;
%       each receives the samples of a pair of classes from the original
%       training dataset, and votes on the most probable class
%
% Options:
%    - subModel: binary classifier used for classification (default: logistic)
%    - subOptions: options for subModel (default: none)
%
% Author:
%  	 - Bita Nejat (2014)

% Parse options
[subModel, subOptions] = myProcessOptions(options, 'subModel', @ml_binaryclass_logistic, 'subOptions', []);

classes = unique(y);
nClasses = length(classes);

% Set up classifiers
nClassifiers = 1/2*nClasses*(nClasses-1);
subModels = cell(nClassifiers,1);
mapping = zeros(nClassifiers,2);

count = 1;
for a = 1:nClasses-1
    for b = a+1:nClasses
        
        % New training data consisting of only class a and b
        XClass1 = X(y==classes(a) ,:);
        XClass2 = X(y==classes(b) ,:);
        yClass1 = y (y==classes(a));
        yClass2 = y (y==classes(b));
        XNew = [XClass1;XClass2];
        yNew = [yClass1; yClass2];
        
        if isfield(subOptions, 'weights')
            zClass1 = subOptions.weights(y==a);
            zClass2 = subOptions.weights(y==b);
            subOptions.weights = [zClass1;zClass2];
        end
        
        % Turn a to 1 and b to -1
        for m = 1:size(yNew)
            if yNew(m) == classes(a)
                yNew(m) = 1;
            else
                yNew(m) = -1;
            end
        end
        
        % Keep class choices for testing
        mapping(count,1) = classes(a);
        mapping(count,2) = classes(b);
        
        subModels{count} = subModel(XNew,yNew,subOptions);
        count = count + 1;
    end
end

% Model outputs
model.mapping = mapping;
model.nClassifiers = nClassifiers;
model.nClasses = nClasses;
model.classes = classes;
model.subModels = subModels;
model.predict = @predict;
model.name = '1-vs-1 Classification';
end

function [yhat] = predict(model,Xhat)
% Prediction Function

[nTest,nFeatures] = size(Xhat);
nClassifiers = model.nClassifiers;
subModels = model.subModels;
mapping = model.mapping;

yhats = zeros(nTest, nClassifiers);
% Run all classifiers
for c = 1: nClassifiers
    
    % Predictions of subModels
    preds = subModels{c}.predict(subModels{c},Xhat);
    
    % Restore actual classes from training
    a = mapping(c,1);
    b = mapping(c,2);
    yhats(preds == 1, c) = a;
    yhats(preds == -1, c) = b;
end

% Calculate and take votes
yhat = mode(yhats, 2);
end
