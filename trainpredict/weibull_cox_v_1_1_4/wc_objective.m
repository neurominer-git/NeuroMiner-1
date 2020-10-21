function [LL, gradient, H] = wc_objective(p, model)

% [LL, gradient, H] = wc_objective(p, model)
%
% Computes the log likelihood of a Weibull-Cox model. Optionally returns
% the gradient and hessian with respect to parameters.
%
% Inputs:
%   p = q+2 by 1 vector where the first q elements correspond to the beta
%       regression coefficients, q+1 = rho, and q+2 = nu.
%   model = model structure
%
% Output:
%   LL = scalar, log likelihood
%   gradient = q+2 by 1 vector where element i gives the partial derivative
%       of the log likelihood with respect to parameter i
%   H = q+2 by q+2 matrix of second order partial derivatives. Element i,j
%       is the second order partial derivaitve of the log likelihood with
%       repect to parameters i and j. (Symmetrix matrix)
%
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

% Extract parameters
beta = p(1:model.q);
beta = beta(:); % ensure column vector
% rho = p(model.q+1);
% nu = p(model.q+2);
if (p(model.q+1) < 700) 
    rho = log(1 + exp(p(model.q+1))); % Ensures non-negativity
    drho = exp(p(model.q+1))/(1 + exp(p(model.q+1)));
else
    rho = p(model.q+1);
    drho = 1;
end
if (p(model.q+2) < 700)
    nu = log(1 + exp(p(model.q+2)));
    dnu = exp(p(model.q+2))/(1 + exp(p(model.q+2)));
else
    nu = p(model.q+2);
    dnu = 1;
end

BX1 = beta'*model.X1';
expBX = exp(beta'*model.X');

% Non-censored individuals
log_lambda = log(nu) - log(rho) + ( (nu-1)*log(model.t1) ) - ( (nu-1)*log(rho) );
L1 = sum(log_lambda) + sum(BX1);

% Censored and non-censored individuals
Lambda = (model.t/rho).^nu;
LexpBX = Lambda.*expBX';
sumLexpBX = sum(LexpBX);
L0 = -sum(LexpBX);

% Prior terms
%Lprior = -0.5*sum(beta'*beta);
Lprior = 0;

LL = -(1/model.N) * ( L0 + L1 + Lprior);

%% Compute gradient
if (nargout > 1)
    
    gradient = zeros(1,length(p));
    % beta
    sumXLexpBX = sum(bsxfun(@times,model.X,LexpBX),1);
    gradient(1,1:model.q) = sum(model.X1,1)...   %L1
        - sumXLexpBX;... %L0
        %- beta'; %Lprior
    
    % rho
    gradrho = (model.N1 * (-nu/rho))... %L1
        + ( (nu/rho)*sumLexpBX );   %L0
    gradient(model.q+1) = drho * gradrho; %due to transformation
    
    % nu
    nuLexpBX = (log(model.t) - log(rho)).*LexpBX;
    gradnu = (model.N1/nu) + sum(log(model.t1/rho))...    %L1
        - sum(nuLexpBX);
    gradient(model.q+2) = dnu * gradnu; %due to transformation
    
    gradient = -(1/model.N) * gradient';
end

%% Compute Hessian
if (nargout > 2)
    
    H = zeros(model.q+2,model.q+2);
    XLexpBX = model.X.*LexpBX;
    
    % beta
    for r = 1:model.q
        H(r,1:model.q) = -sum(XLexpBX.*model.X(:,r));
    end
    H(1:model.q,1:model.q) = H(1:model.q,1:model.q);% - eye(model.q); % prior
    
    % rho
    H(model.q+1,model.q+1) = model.N1*(nu/rho^2) - ( (nu*(nu+1)/rho^2)*sumLexpBX );
    drho2 = drho/(1 + exp(p(model.q+1)));    
    H(model.q+1,model.q+1) = (H(model.q+1,model.q+1)*(drho^2))+...
        (gradrho * drho2);

    % nu
    H(model.q+2,model.q+2) = -model.N1*(1/nu^2)...
        - sum( ( (log(model.t) - log(rho)).^2 ).*LexpBX);
    dnu2 = dnu/(1 + exp(p(model.q+2)));    
    H(model.q+2,model.q+2) = (H(model.q+2,model.q+2)*(dnu^2))+...
        (gradnu * dnu2);

    % rho mu cross term
    H(model.q+1,model.q+2) = -(model.N1/rho)...
        +sum( ((nu/rho)*( log(model.t) - log(rho) ) .* LexpBX) + ( (1/rho)*LexpBX ) );
    H(model.q+1,model.q+2) = (drho * dnu * H(model.q+1,model.q+2));
    H(model.q+2,model.q+1) = H(model.q+1,model.q+2);    % These derivatives are equal
    
    % rho beta cross terms
    H(model.q+1,1:model.q) = (nu/rho)*sumXLexpBX;
    H(model.q+1,1:model.q) = (drho * H(model.q+1,1:model.q));
    H(1:model.q,model.q+1) = H(model.q+1,1:model.q);
    
    % nu beta cross terms
    H(model.q+2,1:model.q) = -sum(model.X.*nuLexpBX);
    H(model.q+2,1:model.q) = (dnu * H(model.q+2,1:model.q));
    H(1:model.q,model.q+2) = H(model.q+2,1:model.q);
    
    H = -(1/model.N) * H;
end

