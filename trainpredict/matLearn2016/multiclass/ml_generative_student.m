function [model] = ml_generative_student(X,~,~)
% ml_multiclass_generativeStudent(X,y,options)
%
% Description:
%	 - Computes the maximum likelihood Student's t fit for the data. Targets
%       y are expected to be empty and are ignored. The prediction function
%       returns an [nTest, 1] vector with the likelihood of each data point
%       in the test set
% 
% Options:
%    - None
%
% Authors:
% 	 - Jason Hartford (2014)
%    - Mark Schmidt

[nTrain, nFeatures] = size(X);

% Optimization parameters
optimOptions.Display = 0;
optimOptions.useMex = 0;

% Calculate mean, sigma and degrees of freedom
vars = minFunc(@studentT,[zeros(nFeatures,1);reshape(eye(nFeatures),nFeatures*nFeatures,1);3],optimOptions,X);
mu = vars(1:nFeatures);
sigma = reshape(vars(nFeatures+1:end-1),nFeatures,nFeatures) + 1e-8*eye(nFeatures);
dof = vars(end);

% Model outputs
model.name = 'Generative Student''s t Model';
model.mu = mu;
model.sigma = sigma;
model.dof = dof;
model.predict = @predict;
end

function [lik] = predict(model,Xhat)
% Predictio function
[nTest, nFeatures] = size(Xhat);
mu = model.mu;
sigma = model.sigma;
dof = model.dof;
Xnorm = bsxfun(@minus, Xhat, mu');

% Form distribution
lik = mvtpdf(Xnorm, reshape(sigma,nFeatures,nFeatures), dof);
end

function [f,g] = studentT(vars,X)
% Student-t loss functon for best fit of mu, sigma and dof

[nTrain,nFeatures] = size(X);
f = 0;
% Derivative wrt mu
g_m = zeros(nFeatures,1);
% Derivative wrt sigma
g_s = zeros(nFeatures);
% Derivative wrt dof
g_d = 0;

mu = vars(1:nFeatures);
sigma = vars(nFeatures+1:end-1);
sigma = reshape(sigma,nFeatures,nFeatures);
dof = vars(end);

[R,err]=chol(sigma);
if err == 0 &&  all(dof > 0)
    sigmaInv = sigma^-1;
    for i = 1:nTrain
        tmp = 1 + (1/dof)*(X(i,:)'-mu)'*sigmaInv*(X(i,:)'-mu);
        f = f + ((nFeatures+dof)/2)*log(tmp);
        
        g_m = g_m - ((nFeatures+dof)/(dof*tmp))*sigmaInv*(X(i,:)'-mu);
        g_s = g_s - ((nFeatures+dof)/(2*dof*tmp))*sigmaInv*(X(i,:)'-mu)*(X(i,:)'-mu)'*sigmaInv;
        g_d = g_d - (nFeatures/(2*tmp*dof^2))*(X(i,:)'-mu)'*sigmaInv*(X(i,:)'-mu);
        g_d = g_d - (dof/(2*tmp*dof^2))*(X(i,:)'-mu)'*sigmaInv*(X(i,:)'-mu);
        g_d = g_d + (1/2)*log(tmp);
    end
else
    f = inf;
    g_m = g_m(:);
    g_s = g_s(:);
    g_d = g_d(:);
    g = [g_m;reshape(g_s,nFeatures*nFeatures,1);g_d];
    return;
end

% Take into account logZ
logSqrtDetSigma = sum(log(diag(R)));
logZ = gammaln((dof+nFeatures)/2) - (nFeatures/2)*log(pi) - logSqrtDetSigma - gammaln(dof/2) - (nFeatures/2)*log(dof);
f = f - nTrain*logZ;
g_s = g_s + (nTrain/2)*sigmaInv;
g_s = (g_s+g_s')/2;
g_d = g_d - (nTrain/2)*psi((dof+nFeatures)/2) + (nTrain/2)*psi(dof/2) + nTrain*(nFeatures/(2*dof));
g = [g_m;reshape(g_s,nFeatures*nFeatures,1);g_d];
end
