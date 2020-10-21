function [model] = ml_binaryclass_multiclass(X,y,options)
% ml_binaryclass_multiclass(X,y,options)
%
% Description:
%	 - Uses a multiclass classifier to treat the special case of binary
%      labels
%
% Options:
%    - subModel: multiclass classifier (default: logistic)
%    - subOptions: options for subModel (default: none)
%    - Bagging: whether randomization with replacement is used in choosing
%       data for subModel (default: 0)
%
% Authors:
% 	 - Radhika Nangia (2014)
  
[nTrain,nFeatures] = size(X);

% Input options
[subModel,subOptions,bagging] = myProcessOptions(options,'subModel', ...
                    @ml_multiclass_logistic,'subOptions',[],'bagging',0);

% Class -1 becomes class 2
y(y == -1) = 2;    

% Randomly reorder and choose data
if bagging
    for i = 1:nTrain
        k = randi(nTrain);
        Xbag(i,:) = X(k,:);
        ybag(i,:) = y(k,:);
    end
    X = Xbag;
    y = ybag;
end

% Model outputs
model.submodel = subModel(X,y,subOptions);
model.name = ['Bin. Class., ', model.submodel.name];
model.predict = @predict;
end

function [yhat] = predict(model,Xhat)
% Predict funtion
yhat = model.submodel.predict(model.submodel, Xhat);

% Change back to class -1
yhat(yhat == 2) = -1;
end