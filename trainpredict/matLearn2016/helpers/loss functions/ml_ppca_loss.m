function [f] = ml_factor_analysis_loss(deltaPhi,nVars,nFactors,S)

delta = reshape(deltaPhi(1:nVars*nFactors),nVars,nFactors);
phi = deltaPhi(end);

sigma = delta*delta' + phi*eye(nVars);

f = logdet(sigma,inf) + trace(sigma\S);