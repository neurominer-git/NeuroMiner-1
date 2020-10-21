function [model] = ml_multiclass_ECOC(X,y,options)
% ml_multiclass_ECOC(X,y,options)
%
% Description:
%    - Relabel classes into binary forms, and use combination of
%       binary classification submodels to predict multiclass labels
%
% Options:
%    - codeDesign: The method of designing coding matrix (default = 'onevsone')
%       'ovo' - one-vs-one
%       'ova' -  one-vs-all
%       'exh' - exhaustive
%       'rnd' - random
%       (default: 'exh')
%    - decodeDesign: The method of decoding the results (default = 'Hamming')
%       'hm' - Hamming
%       'euc' - Euclidean
%       (default: 'hm')
%    - subModel: binary classifier used (default: logistic)
%    - subOptions: options for subModel (default: none)
%
% Authors:
% 	- Saeed Karimifard (2014)

[nTrain, nFeatures] = size(X);

% Load options
[subModel,subOptions,codeDesign,...
 decodeDesign] = myProcessOptions(options,...
                                  'subModel',@ml_binaryclass_logistic,...
                                  'subOptions',[],...
                                  'codeDesign', 'exh',...
                                  'decodeDesign', 'hm');

% Number of unique y labels
classes = unique(y);
nClasses = length(classes);

% Create coding matrix mapping classes to nClassifiers of binary labels
if strcmp(codeDesign, 'rnd')
    nClassifiers = nClasses;
    codingMatrix = (randi([0 1],nClasses,nClassifiers)*2 - 1);
elseif strcmp(codeDesign, 'ovo')
    nClassifiers = (nClasses * (nClasses - 1)) / 2;
    codingMatrix = zeros(nClasses,nClassifiers);
    indx = 1;
    for a = 1:(nClasses - 1)
        for b = (a+1):nClasses
            codingMatrix(a,indx) = 1;
            codingMatrix(b,indx) = -1;
            indx = indx + 1;
        end
    end
elseif strcmp(codeDesign, 'ova')
    nClassifiers = nClasses;
    codingMatrix = -1*ones(nClasses, nClasses) + 2*eye(nClasses);
elseif strcmp(codeDesign, 'exh')
    nClassifiers = (2^(nClasses-1)) - 1;
    codingMatrix = zeros(nClasses,nClassifiers);
    if (nClasses < 3)
        disp('ERROR: Number of classes is not sufficient.');
    else
        codingMatrix(1,:) = 1;
        for i = 2:nClasses
            width = 2^(nClasses-i);
            for k = 1:nClassifiers
                if ( mod(floor((k-1)/width),2) == 0)
                    codingMatrix(i,k) = -1;
                else
                    codingMatrix(i,k) = 1;
                end
            end
        end
    end
end

% Build binary classifier models
for k = 1:nClassifiers
    % Create new training X and y data
    data = zeros(nTrain,1);
    for c =  find(codingMatrix(:,k) == 1)';
        data(y == classes(c)) = 1;
    end
    for c = find(codingMatrix(:,k) == -1)';
        data(y == classes(c)) = -1;
    end
    
    indClass1 = find(data == 1);
    indClass2 = find(data == -1);
    
    XClass1 = X(indClass1, :);
    XClass2 = X(indClass2, :);
    yClass1 = ones(length(indClass1),1);
    yClass2 = -1*ones(length(indClass2),1);
    XNew = [XClass1; XClass2];
    yNew = [yClass1; yClass2];
    
    % Create subModels
    subModels{k} = subModel(XNew,yNew,subOptions);
end

% Model outputs
model.codingMatrix = codingMatrix;
model.nClassifiers = nClassifiers;
model.nClasses = nClasses;
model.classes = classes;
model.predict = @predict;
model.decodeDesign = decodeDesign;
model.subModels = subModels;
model.name = 'Classification using Error-Correcting Output Codes';
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);
nClasses = length(model.classes);

subModels = model.subModels;
codingMatrix = model.codingMatrix;
nClassifiers = model.nClassifiers;

% Predict binary labels using subModels
group = zeros(nTest, nClassifiers);
for c = 1:nClassifiers
    group(:,c) = subModels{c}.predict(subModels{c},Xhat);
end

% Decode actual label by comparing binary predicted labels to codingMatrix
if strcmp(model.decodeDesign, 'hm')
    for i = 1:nTest
        for j = 1:nClasses
            curPred = group(i,:);
            curPred(codingMatrix(j,:) == 0) = [];
            codingVec = codingMatrix(j,:);
            codingVec(codingMatrix(j,:) == 0)  = [];
            dist(i,j) = pdist2(curPred,codingVec,'hamming');
        end
    end
elseif strcmp(model.decodeDesign, 'euc')
    for i = 1:nTest
        for j = 1:nClasses
            curPred = group(i,:);
            curPred(codingMatrix(j,:) == 0) = [];
            codingVec = codingMatrix(j,:);
            codingVec(codingMatrix(j,:) == 0)  = [];
            dist(i,j) = pdist2(curPred,codingVec,'euclidean');
        end
    end
end

% Determine closest actual label
[~, indx] = min(dist, [], 2);
yhat = model.classes(indx);
end
