function r = avrRank(x)
%AVRRANK Returns the sample ranks of the values in a vector. In case of ties
% (i.e., equal values) and missing values it computes their average rank.
%  
% INPUT:
%  x - any row or column vector, may be integer, character
%      or floating point. 
%
% OUTPUT:
%  r - vector of corresponding ranks for each element in x tied elements
%      are given an average rank for all tied elements.
%      
% SEE ALSO:
%  * rankt - package at 
%          http://www.mathworks.com/matlabcentral/fileexchange/27701-rankt
%  * FractionalRankings from Rankings toolbox
%          http://www.mathworks.com/matlabcentral/fileexchange/19496 
%  * nm_tiedrank from MATLAB Statistics Toolbox
%  * rank function in EXCEL
%  * rank function in R language
%
% EXAMPLE:
% x = round(rand(1,10)*5);
% r = avrrank(x);
% disp([x;r])
%     3	  5    2	 1	 5	 0	 3	 2	 5	 0
%    6.5 9.0	4.5	3.0	9.0	1.5	6.5	4.5	9.0	1.5	
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
% Code covered by BSD License

%% Sort and create ranks array
[sx, sortidx] = sort(x(:));
nNaN  = sum(isnan(x));                 % number of NAN's in the array
nx    = numel(x) - nNaN;               % length of x without NAN's
ranks = [1:nx NaN(1,nNaN)]';           % push NAN's to the end

%% Adjust for ties by averaging their ranks
pos = 1;
while pos<nx
  while (pos<nx && sx(pos)~=sx(pos+1)), pos=pos+1; end 
  tieStart = pos;   % search until next tie is found
  while (pos<nx && sx(pos)==sx(pos+1)), pos=pos+1; end
  tieEnd   = pos;   % end of sequence of ties is found
  idx = tieStart : tieEnd;
  ranks(idx) = mean(ranks(idx)); % replace ranks with average ranks
end

%% Prepare output
r(sortidx) = ranks;
r = reshape(r,size(x));                % reshape to match input

