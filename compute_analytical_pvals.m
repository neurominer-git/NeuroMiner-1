function p_map = compute_analytical_pvals(labels, X, w_star)
% p_map = compute_analytical_pvals(labels, X, w_star)
% Compute analytical multivariate P values using the formulation of Gaonkar
% et al. "Interpreting support vector machine models for multivariate group
% wise analysis in neuroimaging", Medical Image Analysis, July 2015.
% =========================================================================
% Nikolaos Koutsouleris, 12/2018

p           = sum(labels==1)/max(size(labels));
r           = size(X,1);
J           = ones(r,1);
K           = X*X';
Z           = inv(K)+(inv(K)*J*inv(-J'*inv(K)*J)*J'*inv(K));
C           = X'*Z;
mean_w      = sum(C,2)*(2*p-1);
SD          = sqrt(sum(C.^2,2)*(4*p-4*p^2)); 
D           = sqrt(sum(SD.*SD+mean_w.*mean_w));
mean_normw  = mean_w/D;
SD_normw    = SD/D;
normw_star  = w_star./norm(w_star);
SD_normw(SD_normw==0)=realmin;
p_map       = 2 * normcdf( -abs( normw_star - mean_normw), zeros( size( w_star ) ), SD_normw);


         