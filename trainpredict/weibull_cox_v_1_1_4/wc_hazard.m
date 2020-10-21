function h = wc_hazard(t,xstar,model)

% h = wc_hazard(t,model)
% 
% Calculates the hazard rate at time t for an individual with
% covariates xstar using a trained model.
% 
% Inputs:
%     t = d-dimensional vector of times at which to evaluate the survival
%         fucntion.
%     xstar = q-dimensional vector of covariates.
%     model = model structure generated using wc_train.
%     
% Output:
%     h = d-dimensional vector of hazard rate values corresponding to
%         the times in vector t.
%
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

xstar = xstar(:);
h = (model.nu/model.rho) * ((t/model.rho).^(model.nu-1)) * exp(model.beta'*xstar);