function X=unmatrizicing(X,n,D)

% Inverse operation of matrizicing, i.e. reconstructs the multi-way array
% from it's matriziced version.
%
% Written by Morten Mørup
%
% Usage:
%   X=unmatrizicing(X,n,D)
%
% Input:
%   D is vector containing original dimensions
%   n is the dimension along which the matrizicing was originally performed
%   X is a matrix
%
% Output:
%   X (multi-way array)
%
% Copyright (C) Morten Mørup and Technical University of Denmark, 
% September 2006
%                                          
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if n==1
    perm=[1:length(D)];
else
    perm=[2:n 1 n+1:length(D)];
end
  
X=permute(reshape(X,D([n 1:n-1 n+1:length(D)])),perm);
