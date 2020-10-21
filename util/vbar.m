function v = vbar(x,y1,y2)
% vbar(x,y1,y2) - creates vertical bar vectors for plotting
%
% If the input arguments are vectors, then vbar returns a complex
% vector v that displays a set of vertical bars (all in a single color)
% when plotted with plt or plot. If the input arguments are matrices
% with n columns, then vbar returns a complex matrix v with n columns,
% so that n sets of vertical bars are plotted in n different colors
% (one color for each column of the inputs).
%
%   if x,y1,y2 are column or row vectors
%     - returns a complex column vector (3 times as long as x)
%
%   if x,y1,y2 are matrices
%     - returns a complex matrix with the same number of columns
%       as x and three times as many rows
%
%   x,y1,y2 all must be the same size.
%   (Although the special case where y1 is scalar is always allowed.)
%
% ----- Author: ----- Paul Mennen
% ----- Email:  ----- paul@mennen.org


[n,nc] = size(x);

if n==1  % special case if one row
  x = x(:); y1 = y1(:); y2 = y2(:); % turn rows into columns
  n = nc;  nc = 1;
end;

if size(y1,1)==1        y1 = repmat(y1,n,1);  end; % expand y1 as needed
if nc>1 & size(y1,2)==1 y1 = repmat(y1,1,nc); end; % expand y1 as needed

v = complex([x; x; x], [y1; y2; x+NaN]);
e = reshape(1:3*n,3,n)';
v(e(:),:) = v;

% Note: The above 3 lines could be replaced by EITHER of the following
% single lines, However in both cases the resulting output is much
% longer than necessary
%       v = quiv(x,y1,0,y2-y1,0);
%       v = ebar(x,y1,0,y2-y1,0);

