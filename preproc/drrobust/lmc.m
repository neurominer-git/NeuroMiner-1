function [res] = lmc(x)

%LMC calculates the left medcouple, a robust measure of
%left tail weight
%
% The left medcouple is described in:
%    Brys, G., Hubert, M. and Struyf, A. (2006),
%    "Robust Measures of Tail Weight",
%    Computational Statistics and Data Analysis,
%    50 (No 3), 733-759. 
%
% For the up-to-date reference, please consult the website:
%    wis.kuleuven.be/stat/robust.html
%
% Required input arguments:
%    x : Data matrix (rows=observations, columns=variables)
%
% I/O:
%    result=lmc(x);
%   
% Example:
%    result = lmc([chi2rnd(5,1000,1) trnd(3,1000,1)]);
%
% The output of LMC is a vector containing the left medcouple
%     for each column of the data matrix x
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Guy Brys
% Last Update: 17/03/2006

if (nargin<1)
    error('No input arguments')
end
if (size(x,1)==1)
    x = x';
end
for (i=1:size(x,2))
    res(i) = -mc(x(x(:,i)<=prctile(x(:,i),50),i));
end
