function [f] = ml_factor_analysis_loss(deltaPhi,nVars,nFactors,S)

delta = reshape(deltaPhi(1:nVars*nFactors),nVars,nFactors);
phi = deltaPhi(nVars*nFactors+1:end);

sigma = delta*delta' + diag(phi);

f = logdet(sigma,inf) + trace(sigma\S);