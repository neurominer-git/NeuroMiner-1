function r=fisherinv(Z)
% function r=fisherinv(Z)
% takes the inverse of the fisher-Z transformation of a correlation coefficients
% A=.5*log((1+A)./(1-A));
%
% see also fisherz
r=(exp(2.*Z)-1)./(exp(2.*Z)+1);
