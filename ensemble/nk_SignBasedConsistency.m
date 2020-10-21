function I = nk_SignBasedConsistency(E)
% Compute modified version of sign-based consistency criterion using a binary
% classifier ensemble. See paper by Vanessa Gomez-Verdejo et al.,
% Neuroinformatics, 2019, 17:593-609. We additionally downweight the
% consistency vector I by the number of nonfinite values in the ensemble matrix.
% I = 2 * abs( nanmean(E>0,2) - 0.5) .* (1 - sum(isnan(E),2) / size( E, 2));
Rnan = 1-sum(isnan(E),2) / size( E, 2 );
Ip = nanmean(E>0,2).*Rnan;
In = nanmean(E<0,2).*Rnan;
I = abs(Ip-In);