function A=fisherz(A)
% function A=fisherz(A)
% takes the fisher-Z transformation of correlation coefficients
% A=.5*log((1+A)./(1-A));
%
% see also fisherinv
A=.5*log((1+A)./(1-A));


