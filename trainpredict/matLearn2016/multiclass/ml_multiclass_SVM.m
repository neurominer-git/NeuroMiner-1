function [model] = ml_multiclass_SVM(X,y,options)
% ml_multiclass_svm(X,y,options)
%
% Description:
%	 - Fits a linear classifier by maximizing the margin using SVM.
%
% Options:
%    - addBias: adds a bias variable (default: 0)
%    - lambdaL2: strenght of L2-regularization parameter (default: 0)
%    - slack: one of
%       'nk' - NK-slack
%       'n' - N-slack
%       (default: 'nk')

[nTrain,nFeatures] = size(X);

% Process options
[addBias,lambdaL2, slack] = myProcessOptions(options,'addBias',1,'lambdaL2',0, 'slack', 'nk');

classes = sort(unique(y));
nClasses = length(classes);

if addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Inputs for linear/quadratic programming
wActual = zeros(nTrain,nClasses*nFeatures);
for i = 1:nTrain
    offset = (find(classes == y(i)) - 1)' * nFeatures;
    wActual(i,offset+1:offset+nFeatures) = 1;
end

wMatrix = (kron(eye(nClasses),ones(nTrain, nFeatures))-repmat(wActual, nClasses, 1));
ind = all(wMatrix == 0, 2);
wMatrix(ind,:) = [];

XMatrix = repmat(X, nClasses, nClasses);
XMatrix(ind,:) = [];

if strcmp(slack, 'n')
    vMatrix = repmat(-eye(nTrain),nClasses,1);
    vMatrix(ind,:) = [];
    A = [XMatrix.*wMatrix,vMatrix; zeros(nTrain,nClasses*nFeatures),-eye(nTrain)];
    b = [-ones((nClasses-1)*nTrain,1); zeros(nTrain,1)];
    f = [zeros(nFeatures*nClasses,1); ones(nTrain, 1)];
    
    options_ = [];
    options_.Display = 'off';
    if lambdaL2 == 0
        x = linprog(f,A,b,[],[],[],[],[],options_);
    else
        % Quadratic programming with L2 regularization
        H = [lambdaL2*eye(nFeatures*nClasses), zeros(nFeatures*nClasses,nTrain); zeros(nTrain, nFeatures*nClasses), zeros(nTrain, nTrain)];
        x = quadprog(H,f,A,b,[],[],[],[],[],options_);
    end
elseif strcmp(slack, 'nk')
    vMatrix = -eye(nClasses*nTrain);
    vMatrix(ind,:) = [];
    
    A = [XMatrix.*wMatrix,vMatrix; zeros(nClasses*nTrain, nClasses*nFeatures), eye(nClasses*nTrain)];
    b = [-ones((nClasses - 1)*nTrain,1); zeros(nClasses*nTrain,1)];
    vMatrix2 = ones(nClasses*nTrain, 1);
    vMatrix2(ind) = 0;
    f = [zeros(nFeatures*nClasses,1); vMatrix2];

    options_ = [];
    options_.Display = 'off';
    if lambdaL2 == 0
        % Linear programming without L2 regularization
        x = linprog(f,A,b,[],[],[],[],[],options_);
    else
        % Quadratic programming with L2 regularization
        H = [lambdaL2*eye(nFeatures*nClasses), zeros(nFeatures*nClasses,nTrain*nClasses); zeros(nTrain, nFeatures*nClasses), zeros(nTrain, nTrain*nClasses)];
        x = quadprog(H,f,A,b,[],[],[],[],[],options_);
    end
end
w = reshape(x(1:nFeatures*nClasses),nFeatures,nClasses);

% Model outputs
model.w = w;
model.addBias = addBias;
model.classes = classes;
if strcmp(slack, 'n')
    model.name = 'N-Slack SVM Classification';
elseif strcmp(slack, 'nk')
    model.name = 'NK-Slack SVM Classification';
end
model.predict = @predict;
end

function [yhat] = predict(model, Xhat)
% Model outputs
[nTest, nFeatures] = size(Xhat);

if model.addBias
    Xhat = [ones(nTest,1) Xhat];
end

% Predict labels (1..nClasses) for xhat
prob = Xhat * model.w;
[~, indx] = max(prob, [], 2);
yhat = model.classes(indx);
end
