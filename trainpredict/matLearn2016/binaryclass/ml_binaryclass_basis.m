function [model] = ml_binaryclass_basis(X,y,options)
% ml_binaryclass_basis(X,y,options)
%
% Description:
%	 - Binary classification after a change of basis
%
% Options:
%	 - subModel: the model used for classification (default: brokenStump)
%    - subOptions: options for subModel (default: none)
%    - basisFunc: function that will perform the basis transformation on the variable X
%       (default: RBF)
%    - basisOptions: options for basisFunc (default: sigma = 1)
%
% Authors:
%    - Daniel Fugere (2014)

% Set default options if not specified by input
basisOptDefault = [];
basisOptDefault.sigma = 1;
[subModel, subOptions, basisFunc, basisOptions] = myProcessOptions(options,'subModel',@ml_binaryclass_brokenStump, 'subOptions', [], 'basisFunc',@ml_kernel_rbf, 'basisOptions', basisOptDefault);

% Calculate new X values after basis change
model.base = X;
[X, name] = basisFunc(X,X,basisOptions);
newX  = X;

% Create model using new X values and specified sub-model
baseModel = subModel(newX,y,subOptions);

% Model outputs
model.name = ['Binary Classification under Basis Change, ', name];
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
