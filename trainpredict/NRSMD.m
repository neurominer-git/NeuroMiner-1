% =========================================================================
% FORMAT param = NMRSD(expected, predicted)
% =========================================================================
% Compute Normalized Root of Mean Squared Deviation (NRMSD) of regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function param = NRSMD(expected, predicted)

sumerr = sum((expected-predicted).^2);
sumexp = sum(expected.^2);
param = sqrt(sumerr/sumexp)*100;

end