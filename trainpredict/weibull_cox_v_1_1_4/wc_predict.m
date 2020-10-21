function [mean_prediction, var_prediction] = wc_predict(xstar,model)

% [mean_prediction, var_prediction] = wc_predict(xstar,model)
%
% Computes the expected event time and variance for an indivdual with
% covariates xstar by numerically computing the mean and variance of the
% associcated event time probability density.
%
% Inputs:
%     xstar = q-dimensional vector of covariates.
%     model = model structure generated using wc_train.
%
% Output:
%     mean_prediction = predicted event time
%     var_prediction = variance of predicted event time
% 
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

THRESHOLD = 1e-12;
LB = 0; % lower integration bound

xstar = xstar(:);   % ensure column vector

if (length(xstar) ~= model.q)
    error('Dimension of xstar must be the same as training data')
end

% Upper Bound
% This begins at the mode of the distribution and takes jumps given by
% 10*variance(event times) until the integrand is sufficiently small.

% Initialises upper bound to the maximum of the event time density. If the
% maximum doesn't exist then it is set to 1.
if ((model.rho^(model.nu))*(1/model.nu)*(model.nu-1)*exp(-model.beta'*xstar) <= 0)
    ub = 1;
else
    ub = ((model.rho^(model.nu))*(1/model.nu)*(model.nu-1)*exp(-model.beta'*xstar))^(1/model.nu);
end
    
ub = ub + (10*model.time_var);
pub = wc_mean_crude(ub,xstar,model);
while (pub > THRESHOLD)
    ub = ub + model.time_var;
    pub = wc_mean_crude(ub,xstar,model);
end

mean_prediction = integral(@(r)wc_mean_crude(r,xstar,model),LB,ub,'RelTol',1e-8);

if (nargout > 1)
    
    if ((model.rho^(model.nu))*(1/model.nu)*(model.nu-1)*exp(-model.beta'*xstar) <= 0)
        ub = 1;
    else
        ub = ((model.rho^(model.nu))*(1/model.nu)*(model.nu-1)*exp(-model.beta'*xstar))^(1/model.nu);
    end
    
    ub = ub + (10*model.time_var);
    pub = wc_mean_crude(ub,xstar,model);
    while (pub > THRESHOLD)
        ub = ub + model.time_var;
        pub = wc_var_crude(ub,xstar,model);
    end
    
    var_prediction = integral(@(r)wc_var_crude(r,xstar,model),LB,ub,'RelTol',1e-8);
    
    var_prediction = var_prediction - (mean_prediction^2);
end