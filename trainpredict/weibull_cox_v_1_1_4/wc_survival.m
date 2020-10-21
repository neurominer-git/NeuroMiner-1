function S = wc_survival(t,xstar,model)

% S = wc_survival(t,model)
% 
% Calculates the survival probability at time t for an individual with
% covariates xstar using a trained model.
% 
% Inputs:
%     t = d-dimensional vector of times at which to evaluate the survival
%         fucntion.
%     xstar = q-dimensional vector of covariates.
%     model = model structure generated using wc_train.
%     
% Output:
%     S = d-dimensional vector of survival probabilities corresponding to
%         the times in vector t.
%
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

xstar = xstar(:);
S = exp( - ((t/model.rho).^model.nu) * exp(model.beta'*xstar) );