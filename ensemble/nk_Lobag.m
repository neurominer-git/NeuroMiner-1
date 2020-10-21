% =========================================================================
% FORMAT [mED, ED] = nk_Lobag(E, L)
% =========================================================================
% This function performs the bias / variance decomposition of the error 
% proposed by Valentini & Dietterich (2004, Journal of Machine Learning 
% Research (5), 725-775) for a given ensemble of base learners' hypotheses 
% (E) with respect to the supervision information of given test samples (L).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2011
function [mED, ED] = nk_Lobag(E, L)

nE = size(E,1);
E = sign(E);
% Bias computation for matrix mode
sum_x       = sum(E,2); 
ind_x       = sum_x > 0;
g_x         = zeros(nE,1);
g_x(ind_x)  = 1; 
g_x(~ind_x) = -1;

% Biased results
bias = g_x ~= L;

% Biased variance 
vb = var(E(bias,:),1,2); 
if isempty(vb), vb = 0; end

% Unbiased variance
vu = var(E(~bias,:),1,2); 
if isempty(vu), vu = 0; end
    
% Lobag computation
ED = mean(bias) + mean(vu) - mean(vb);

mED = mean(ED);

end