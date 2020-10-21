function [model] = ml_regression_basis(X,y,options)
% ml_regression_basis(X,y,options)
%
% Description:
%    - Regression under a change of basis
%
% Options:
%	 - subModel: The model used for regression (default: mean)
%    - subOptions: Options for submodel (default: none)
%    - basisFunc: Function for performing basis change on X
%       (default: polynomial1D)
%    - basisOptions: Options for basisFunc (default: polynomial degree <= 2)
%
% Authors:
%    - Roee Bar (2014)

% Set default options if not specified by input
basisOptDefault = [];
basisOptDefault.sigma = 1;
[subModel,basisFunc, subOptions, basisOptions] = myProcessOptions(options,'subModel',@ml_regression_mean,'basisFunc',@ml_kernel_rbf, 'subOptions', [], 'basisOptions', basisOptDefault);

% Calculate new X values after basis change
model.base = X;
[X, name] = basisFunc(X,X,basisOptions);
newX  = X;

% Create model using new X values and specified sub-model
baseModel = subModel(newX,y,subOptions);

% Model outputs
model.name = ['Regression under Basis Change with: ', name];
model.newX = newX;
model.basisFunc = basisFunc;
model.basisOptions = basisOptions;
model.baseModel = baseModel;
model.predict = @predict;

end

function [yhat] = predict(model,Xhat)
% Prediction function
% Input basis-changed Xhat into submodel's predict function
newX = model.basisFunc(Xhat,model.base,model.basisOptions);
yhat = model.baseModel.predict(model.baseModel,newX);
end