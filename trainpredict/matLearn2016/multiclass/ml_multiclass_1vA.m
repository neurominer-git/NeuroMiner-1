function [model] = ml_multiclass_1vA(X,y,options)
% ml_multiclass_1vA(X,y,options)
%
% Description:
%    - Classification using one classifier for each class, labelled as 1, against all other classes,
%	labelled as 0. Note that learning on balanced classifiers will reflect unbalanced
%       distributions. Due to the nature of setting 1 for a given class 
%       and 0 to ALL, this trains each binary classifier towards a negatively skewed 
%       label distribution.
%
% Options:  
%	 - subModel: binary classifier where output is a vector of probabilities
%       (default: logistic)
%	 - subOptions:
%
% Author:
% 	-  Xi Laura Cang (2014) 

[nTrain,nFeatures] = size(X);

% Parse Options
[subModel, subOptions] = myProcessOptions(options, 'subModel', @ml_binaryclass_logistic, 'subOptions', []);

% Number of unique classes of y
classes = unique(y);
nClasses = length(classes);

% Set up classifiers
subModels = cell(nClasses);
for c = 1:nClasses
    yNew = zeros(nTrain, 1);
    for n = 1:nTrain
       if y(n) == classes(c)
           % Label all points in class c as 1
           yNew(n) = 1;
       else
           % Label all other points as -1
           yNew(n) = -1;
       end
    end
    subModels{c} = subModel(X, yNew, subOptions);
end    

% Model outputs
model.subModels = subModels;
model.predict = @predict;
model.nClasses = nClasses;
model.classes = classes;
model.name = '1-vs-All Classification';
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);

subModels = model.subModels;
nClasses = model.nClasses;

yhats = zeros(nTest, nClasses);
% Run all classifiers
for c = 1:nClasses
    probs = subModels{c}.predictProb(subModels{c},Xhat);
    yhats(:,c) = probs(:,2);
end

% Pick the most probable class
[~, indx] = max(yhats, [], 2);
yhat = model.classes(indx);
end
