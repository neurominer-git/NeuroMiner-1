function [model] = ml_regression_totalL2(X,y,options)
% ml_regression_totalL2(X,y,options)
%
% Description:
%    - Returns total least squares minimizing L2 in both X and y
%
% Options:
%    - Bias: adds a bias variable (default: 1)

[nTrain,nFeatures] = size(X);

% Set default options if not specified by input
[addBias] = myProcessOptions(options,'addBias',1);

% Add bias variable if required
if addBias
		X = [ones(nTrain,1) X];
		nFeatures = nFeatures + 1;
end

% Singular value decompositition method of calculating weights
[U,D,V] = svd([X y]);
w = -V(1:nFeatures,nFeatures+1)/V(nFeatures+1,nFeatures+1);

% Model outputs
model.w = w;
model.predict = @predict;
model.name = 'Total Least Squares Regression';
model.addBias = addBias;
end

function [yhat] = predict(model,Xhat)
% Prediction function
	[nTest,nFeatures] = size(Xhat);
	if model.addBias
		Xhat = [ones(nTest,1) Xhat];
        nFeatures = nFeatures + 1;
	end
	yhat = Xhat*model.w;
end
