function model =  wc_train(X, t, E, optionsdef)

% model =  wc_train(X, t, E)
% 
% Fits a Cox proportional hazards model.
% 
% Inputs:
%     X = N by q matrix of covariates where each row corresponds to an
%         observation and each column corresponds to a covariate.
%     t = N by 1 vector of event times.
%     E = N by 1 vector of indicator variables where '0' corresponds to 
%         right censoring and '1' corresonds to the primary event.
%     
% Output:
%     model = trained model structure. In particular model.beta are the
%             regression coefficitions and model.beta_std are the
%             correspoinding error bars. Similarly, model.rho and model.nu
%             are the inferred values for rho and nu. The rest of the items
%             are for internal use.
%
% See weibull_cox_doc_v_1_1_2.pdf for more information.
%
% Run wc_example.m for an illustration of how to use the functions.
%
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

JITTER = 1e-2;  % Ensures safe inversion of hessian

% Some error checks

if (~isvector(t))
        error('Event times must be a column vector')
end
if (~isvector(E))
        error('indicator variables must be a column vector')
end
if (~ismatrix(X))
    error('X must be a matrix')
end

% Set options for optimisation
options = optimset;
options.GradObj = 'on';
options.Hessian = 'on';
options.Display = 'off';
options.FunValCheck = 'off';
options.UseParallel = true;
if ~exist('optionsdef','var')
    options.MaxIter = 10000;
    options.MaxFunEvals = 10000;
    options.TolX = 1e-12;
    options.TolFun = 1e-12;
    options.Algorithm = 'trust-region-reflective';
else
    options.MaxIter = optionsdef.MaxIter;
    options.MaxFunEvals = optionsdef.MaxFunEvals;
    options.TolX = optionsdef.TolX;
    options.TolFun = optionsdef.TolFun;
    options.Algorithm = optionsdef.Algorithm;
end
% Create model structre
[model.N, model.q] = size(X);

if (length(t) ~= model.N)
        error('Event time vector must have the same no. of rows as covariate matrix')
end

if (length(E) ~= model.N)
        error('Indicator vector must have the same no. of rows as covariate matrix')
end

ind = find(E==0);
model.X = X;
model.t = t;
model.E = E;
model.X0 = X(ind,:);
model.t0 = t(ind);
model.N0 = length(model.t0);

ind = find(E==1);
model.X1 = X(ind,:);
model.t1 = t(ind);
model.N1 = length(model.t1);
model.time_var = var(model.t);  % used for a rough time scale in predictions

% Optimise parameters
initial_p = ones(model.q+2,1);   % [beta, rho, nu]
[p, model.minLL] = fminunc(@(par)wc_objective(par,model), initial_p, options);
model.beta = p(1:model.q);
if (p(model.q+1) < 700)
    model.rho = log(exp(p(model.q+1)) + 1 );
else
    model.rho = p(model.q+1);
end
if (p(model.q+2) < 700)
    model.nu = log(exp(p(model.q+2)) + 1 );
else
    model.nu = p(model.q+2);
end

% Compute error bars
% [~, ~, H] = wc_objective(p, model);
% H = (model.N*(H + H')/2) + JITTER*eye(model.q+2);
% L = chol(model.N*H);    % Use cholesky decomposition to invert A
% invL = inv(L);
% invH = invL*invL';
% errs = sqrt(diag(invH));
% model.beta_error = errs(1:model.q);
% model.cov_par = invH;