function [n, nn] =nm_nancount(Y)
% function n=nm_nancount(Y)
n=sum(~isnan(Y)==1);
nn=sum(isnan(Y)==1);