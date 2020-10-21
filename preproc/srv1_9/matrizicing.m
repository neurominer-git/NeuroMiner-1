function Y=matrizicing(X,n)

% The matrizicing operation, i.e.
% matrizicing(X,n) turns the multi-way array X into a matrix of dimensions:
% (size(X,n)) x (size(X,1)*...*size(X,n-1)*size(X,n+1)*...*size(X,N))
% where N=ndims(X)
%
% Input:
% X     Multi-way array
% n     Dimension along which matrizicing is performed
%
% Output:
% Y     The corresponding matriziced version of X
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

N=ndims(X);
Y=reshape(permute(X, [n 1:n-1 n+1:N]),size(X,n),prod(size(X))/size(X,n));
