% centroid.m
% Moody T. Chu
%
% Slightly modified by Jason R. Blevins
%
% $Id: centroid.m,v 1.3 2004/07/15 19:21:56 jrblevin Exp $


function [B,V,Z] = centroid(A);

num_iteration = 100;

[n,m] = size(A);

Aold = A;
B = [];
V = [];
Cv = [];
Z = [];

p = rank(A);
%% p = 100;

%% flops(0)

for i = 1:min(p,num_iteration)
   [A_next,z,b,v,npmu] = new_decom(A);
   B = [B,b];
   Z = [Z,z];
   V = [V,v];
   Cv = [Cv,npmu];
   A = A_next;
end
%% disp('Total flops used = '),flops

% Note that Rold = R_next + Z*inv(diag(S))*Z';
%
% Note also that S/n is the approximate singular value

function [A_next,z,b,v,npmu] = new_decom(A);
%
% This code is for Algorithm 8.1
%
% Input:
%
%		A = indexing matrix (n-by-m)
%
% Output:
%
%		npmu = n*(centroid values)
%		v = centroid factors
%		b = loadings
%

[n,m] = size(A);
tol = eps*max([n,m]);

for i = 1:n
   d(i,1) = norm(A(i,:))^2;				% diagonal of AA'
end

% this is the start value

z = ones(n,1);
k = 1;
w = A'*z;

while ~isempty(k)
   
   e = d - sign(z).*(A*w);				% Step 1
   [y,k] = max(e);					% Step 2
   
   if y > 0 & abs(y) > tol    
      z(k) = -z(k);					% Step 3
      if z(k) == 1					% Step 4
         w = w + 2*A(k,:)';
      else
         w = w - 2*A(k,:)';
      end
   else
      k = [];
   end
end

npmu = w'*w;
v = w/sqrt(npmu);					% Centroid factors
b = A*v;						% Loadings
A_next = A - b*v';