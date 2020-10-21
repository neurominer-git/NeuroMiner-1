function [model] = ml_binaryclass_SVM(X,y,options)
% ml_binaryclass_svm(X,y,options)
%
% Description:
%	 - Fits a linear classifier by maximizing the margin using SVM
%
% Options:
%    - addBias: adds a bias variable (default: 0)
%    - lambdaL2: strenght of L2-regularization parameter (default: 0)

[nTrain,nFeatures] = size(X);

% Process options
[addBias,lambdaL2] = myProcessOptions(options,'addBias',1,'lambdaL2',0);

if addBias
   X = [ones(nTrain,1) X];
   nFeatures = nFeatures + 1;
end

% Inputs for linear/quadratic programming
A = [-X.*repmat(y, [1, nFeatures]), -eye(nTrain); zeros(nTrain, nFeatures), -eye(nTrain)];
b = [-ones(nTrain,1); zeros(nTrain,1)];
f = [zeros(nFeatures,1); ones(nTrain, 1)];
options_ = [];
options_.Display = 'off';

if lambdaL2 == 0
    % Linear programming without L2 regularization
    x = linprog(f,A,b,[],[],[],[],[],options_);
else
    % Quadratic programming with L2 regularization
    H = [lambdaL2*eye(nFeatures), zeros(nFeatures,nTrain); zeros(nTrain, nFeatures), zeros(nTrain, nTrain)];
    x = quadprog(H,f,A,b,[],[],[],[],[],options_);
end
try
    w = x(1:nFeatures);
catch
    error('Linprog did not converge! Use some degree of L2 regularization');
end
% Model outputs
model.w = w;
model.addBias = addBias;
model.name = 'SVM Binary Classification';
model.predict = @predict;
end

function p = predict(model, Xhat)
% Prediction function
[nTest, nFeatures] = size(Xhat);

if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end

p.D = Xhat*model.w; p.D(isnan(p.D))=0; 
index =  p.D >= 0;
p.yhat = sign(index - 0.5);

end