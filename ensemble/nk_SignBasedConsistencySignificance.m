function [P, Pfdr, Z, I] = nk_SignBasedConsistencySignificance(E)

% Computed sign-based consistency [0= no consistency, 1=absolute
% consistency]
I = nk_SignBasedConsistency(E);
% Compute Z score
Z = I ./ std(I);
% Compute P value of Z score (right-tailed one-sided test of Z)
P = normcdf(-Z,0,1);
% Correct for multiple comparisons using the FDR-rate
[~,~,~,Pfdr] = fdr_bh(P,0.05);
%log-transform P values
P = -log10(P); Pfdr = -log10(Pfdr);