function q = quiv(A,B,h,in4,in5)
% If the input arguments are vectors, then quiv returns a complex vector
% that displays a vector field when plotted with plt or plot.
% If the input arguments are matrices, then quiv returns a complex matrix
% that displays a series of vector fields when plotted with plt or plot
% (one vector field per column of the input matrices).
%
% There are 8 possible calling sequences for quiv depending on whether the
% input arguments are real or complex and on whether the optional arrow
% head size argument is included. Quiv is smart enough to figure out which
% calling sequence you are using.
%
% Calling sequence     Tail coordinates    Arrow width/length   Arrow head size
% ------------------   -----------------   ------------------   ---------------
% q = quiv(A,B)        [real(A) imag(A)]   [real(B) imag(B)]    0.3
% q = quiv(A,B,h)      [real(A) imag(A)]   [real(B) imag(B)]    h
% q = quiv(x,y,B)      [x y]               [real(B) imag(B)]    0.3
% q = quiv(x,y,B,h)    [x y]               [real(B) imag(B)]    h
% q = quiv(A,u,v)      [real(A) imag(A)]   [u v]                0.3
% q = quiv(A,u,v,h)    [real(A) imag(A)]   [u v]                h
% q = quiv(x,y,u,v)    [x y]               [u v]                0.3
% q = quiv(x,y,u,v,h)  [x y]               [u v]                h
%
% where:
%   A,B are complex vectors or matrices
%   x,y,u,v are real vectors or matrices
%   h is a scalar
%
% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org

hDef = 0.3;  % default arrowhead size (relative to arrow length)
if nargin==2 h=hDef;
else
  switch  nargin + isreal(A) + isreal(B) - 3*isreal(h)
    case  1,                   B = complex(B,h);   h=hDef;  % q = quiv(A,u,v)
    case  2,                   B = complex(B,h);   h=in4;   % q = quiv(A,u,v,h)
    case  3, A = complex(A,B); B = complex(h,in4); h=hDef;  % q = quiv(x,y,u,v)
    case  4, A = complex(A,B); B = complex(h,in4); h=in5;   % q = quiv(x,y,u,v,h)
    case  5, A = complex(A,B); B = h;              h=hDef;  % q = quiv(x,y,B)
    case  6, A = complex(A,B); B = h;              h=in4;   % q = quiv(x,y,B,h)
  end;
end;

% q = quiv(A,B,h) -------------------------------------------------
e = complex(h,h*h);
[nr nc] = size(A);
n = min(nr,nc);
if n==1  A=A(:); B=B(:); end;  % transform row vectors to column vectors
T = A + B;    % T = position of arrow head tip
P = T-B*e;
Q = T-B*e';
N = A + NaN;
% note: we don't really need the special case below
%       for n==1, but it is probably more efficient.
if n==1 q = reshape(transpose([A T N P T Q N]),[],1);             % single column
else    q = reshape(permute(cat(3,A,T,N,P,T,Q,N),[3 1 2]),[],nc); % multiple columns
% Note: the line above is a vectorized form of this code ---
%   q = [];
%   N = A(:,1) + NaN;
%   for k=1:nc
%      v = transpose([A(:,k) T(:,k) N P(:,k) T(:,k) Q(:,k) N]);
%      q = [q v(:)];
%   end;
end;
