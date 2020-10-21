function [model] = ml_multiclass_basis(X, y, options)
% ml_multiclass_basis(X,y,options)
%
% Description:
%	 - Classification with a basis change.
%
% Options:
%    - subModel: model used for classification (Defaulut: GDA)
%    - subOptions: options for subModel (default: none)
%    - basisFunc: function used to transform X (default: RBF)
%    - basisOptions: options for basisFunc (default: sigma = 1);
%
% Authors:
%    - Daniel Fugere (2014)

if nargin < 3
    fprintf('Only %f input arguments found, please provide all 3: X,y and options as shown in the example',nargin);
end

basisOptDefault = [];
basisOptDefault.sigma = 1;
[subModel, subOptions, basisFunc, basisOptions] = myProcessOptions(options,'subModel',@ml_multiclass_generativeGDA, 'subOptions', [], 'basisFunc',@ml_kernel_rbf, 'basisOptions', basisOptDefault);

% Calculate new X values after basis change
model.base = X;
[X, name] = basisFunc(X,X,basisOptions);
newX  = X;

% Create model using new X values and specified sub-model
baseModel = subModel(newX,y,subOptions);

% Model outputs
model.name = ['Classification under Basis Change with: ', name];
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
