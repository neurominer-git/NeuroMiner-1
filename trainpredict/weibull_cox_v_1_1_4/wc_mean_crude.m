function p = wc_mean_crude(r, xstar, model)

% p = wc_mean_crude(s, xstar, model)
% 
% Copyright James Barrett 2014
%
% Version 1.1.4
% Date: 21 July 2014
% Contact: james.j.barrett@kcl.ac.uk

xstar = xstar(:);

lambda = (model.nu/model.rho) * ((r/model.rho).^(model.nu-1));
Lambda = (r/model.rho).^model.nu;

expBX = exp(model.beta'*xstar);

p = r .* lambda .* (expBX * exp(-Lambda*expBX));


