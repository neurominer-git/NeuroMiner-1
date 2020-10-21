function [model] = ml_multiclass_generativeGDA(X,y,options)
% ml_multiclass_generativeGDA(X,y,options)
%
% Description:
%	 - Classification using a Gaussian generative model with shared full covariance matrices to fit each class
%
% Options:
% 	 - subModel: generative model to fit each class (default: Gaussian)
%    - subOptions: options to subModel (default: none)
%
% Authors:
% 	 - Mark Schmidt
%    - Neil Traft (2014)

[nTrain, nFeatures] = size(X);

% Default options
[subModel subOptions] = myProcessOptions(options, 'subModel', @ml_generative_Gaussian, 'subOptions',[]);

% Number of unique classes of y
classes = unique(y);
nClasses = length(classes);

% Fit generative model to each class
for i = 1:nClasses
    X_c = X(y==classes(i),:);
    nTrain_c = size(X_c,1);
    priors(i) = nTrain_c / nTrain;
    subModels{i} = subModel(X_c,y,subOptions);
end

% Model outputs
model.subModels = subModels;
model.priors = priors;
model.predict = @predict;
model.probability = @predictProb;
model.nClasses = nClasses;
model.classes = classes;
model.name = ['Discr. Classification: ', subModels{1}.name];
end

function [phat] = predictProb(model,Xhat)
% Prediction probabibily function
[nTest, nFeatures] = size(Xhat);

% Determine probability of obtaining each class using fitted generative models
phat = zeros(nTest,model.nClasses);
for i = 1:model.nClasses
    phat(:,i) = model.priors(i)*model.subModels{i}.predict(model.subModels{i},Xhat);
end
phat = bsxfun(@rdivide, phat, sum(phat,2));
end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);
phat = predictProb(model,Xhat);

% Determine class of maximum probability
[~, indx] = max(phat,[],2);
yhat = model.classes(indx);
end
