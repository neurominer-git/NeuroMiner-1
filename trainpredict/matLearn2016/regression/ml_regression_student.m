function [model] = ml_regression_student(X,y,options)
% ml_regression_student(X,y,options)
%
% Description:
%	 - implements linear (or polynomial regression) based on
%       student-t loss function
%
% Options: options for the regression model, with the following fields:
%    - weights: vector of weights for training data instances
%       (default: all 1's)
%    - lambdaL2: values >=0 for L2 regularization (default: 0)
%    - addBias: accepts 0 or 1. If 1, adds bias to X (default: 1)
%    - poly: values >=1 which determines the polynomial
%       order for polynomial regression (default: p=1; linear regression)
%
% Outputs:
%    - model.w: nFeature*1 vector of regression coefficients
%    - model.v: degrees of freedom for the student-t distribution
%    - model.sig2inv: equals 1/sigma^2, sensitivity parameter of student-t dist
%
%    The following fields are the same as the fields in options:
%    - model.addBias, model.lambdaL2, model.weights, model.polyand
%
% Author:
%    - Seyed Ali Saberali (2014)

[nTrain,nFeatures] = size(X);

% Check the existing fields in options and set default values for parameters not specified by user
[addBias, lambdaL2, polyOrder, z] = ...
    myProcessOptions(options, 'addBias', 1, ...
    'LambdaL2', 0, 'poly', 1, 'weights', ones(nTrain, 1));

% Generate the design matrix
if polyOrder>1
    X = repmat(X,1,polyOrder+1).^(repmat(0:polyOrder,size(X,1),1));
    nFeatures = nFeatures + polyOrder+1;
elseif addBias
    X = [ones(nTrain,1) X];
    nFeatures = nFeatures + 1;
end

% Optimization options
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Initial value for the sensitivity parameter of student distribution
sig2Inv=1;
% Initial value for the degrees of freedom of student distribution
v=2;
% Initialize the regression coefficient vector
w0=rand(nFeatures,1);

if lambdaL2 == 0
    % Optimize parameters without regularization
    vars = minFunc(@student_tLoss,[w0;sig2Inv;v],optimOptions,X,y,z);
else
    % Optimize parameters with L2 regularization
    vars = minFunc(@student_tLossL2,[w0;sig2Inv;v],optimOptions,X,y,z,lambdaL2);
end

% Model Outputs
model.name = 'Student''s Loss Regression';
model.w = vars(1:end-2);
model.sig2Inv = vars(end-1);
model.v = vars(end);
model.addBias = addBias;
model.poly=polyOrder;
model.predict = @predict;

end

function [yhat] = predict(model,Xhat)
% Prediction function
[nTest,nFeatures] = size(Xhat);
if model.poly>1
    Xhat = repmat(Xhat,1, model.poly+1).^(repmat(0: model.poly,size(Xhat,1),1));
    
elseif model.addBias
    Xhat = [ones(nTest,1) Xhat];
    nFeatures = nFeatures + 1;
end
yhat = Xhat*model.w;
end

function [f,g] = student_tLoss(vars,X,y,z)
% Student-t Loss function without L2 regularization
% Univariate student PDF used is given byRet
% p(y;sig2In,v,u)=v^(v/2)*gamma(v+1/2)/pi^1/2/gamma(v/2)*sqrt(sig2Inv)...
% *(v+(y-u)^2*sig2Inv)^(-v/2-1/2)

% Weights
w=vars(1:end-2);
% Precision parameter of student distribution
sig2Inv=vars(end-1);
% Degrees of freedom for student distribution
v=vars(end);

if (sig2Inv<0)||(v<0)
    f=inf;
    g=zeros(length(vars),1);
else
    C=v/2*log(v)+log(gamma((v+1)/2))-1/2*log(pi)-log(gamma(v/2))+1/2*log(sig2Inv);
    Xw_y=X*w-y;
    V=-(v+1)/2*log(v+Xw_y.^2*sig2Inv);
    f=-C*sum(z)-sum(z.*V);
    
    R=Xw_y*sig2Inv./(v+Xw_y.^2*sig2Inv);
    g_w=(v+1)*X'*(z.*R);
    g_sig2Inv=-sum(z)/2/sig2Inv+(v+1)/2*sum(z.*(Xw_y.^2./(v+sig2Inv*Xw_y.^2)));
    g_v=-sum(z)/2*(log(v)+1+psi((v+1)/2)-psi(v/2))+...
        1/2*sum(z.*(log(v+Xw_y.^2*sig2Inv)+(v+1)./(v+Xw_y.^2*sig2Inv)));
    g=[g_w;g_sig2Inv;g_v];
end
end

function [f,g] = student_tLossL2(vars,X,y,z,lambda)
% Student-t Loss function with L2 regularization
w=vars(1:end-2);
sig2Inv=vars(end-1);
v=vars(end);

if (sig2Inv<0)||(v<0)
    f=inf;
    g=zeros(length(vars),1);
else
    C=v/2*log(v)+log(gamma((v+1)/2))-1/2*log(pi)-log(gamma(v/2))+1/2*log(sig2Inv);
    Xw_y=X*w-y;
    V=-(v+1)/2*log(v+Xw_y.^2*sig2Inv);
    % Only apply regularization to regression coefficients;
    b=[w;0;0];
    f=-C*sum(z)-sum(z.*V)+(lambda/2)*(b'*b);
    
    R=Xw_y*sig2Inv./(v+Xw_y.^2*sig2Inv);
    g_w=(v+1)*X'*(z.*R);
    g_sig2Inv=-sum(z)/2/sig2Inv+(v+1)/2*sum(z.*(Xw_y.^2./(v+sig2Inv*Xw_y.^2)));
    g_v=-sum(z)/2*(log(v)+1+psi((v+1)/2)-psi(v/2))+...
        1/2*sum(z.*(log(v+Xw_y.^2*sig2Inv)+(v+1)./(v+Xw_y.^2*sig2Inv)));
    g=[g_w;g_sig2Inv;g_v] + lambda*b;
end
end
