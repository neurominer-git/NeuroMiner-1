function R=nk_kurtosis(i,DIM)
% KURTOSIS estimates the kurtosis
%
% y = kurtosis(x,DIM)
%   calculates kurtosis of x in dimension DIM
%
% DIM	dimension
%	1: STATS of columns
%	2: STATS of rows
%	default or []: first DIMENSION, with more than 1 element
%
% features:
% - can deal with NaN's (missing values)
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, VAR, STD, VAR, SKEWNESS, MOMENT, STATISTIC, 
%    IMPLICIT_SKIP_NAN
%
% REFERENCE(S):
% http://mathworld.wolfram.com/

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%	Copyright (C) 2000-2003 by Alois Schloegl <a.schloegl@ieee.org>	
%	$Revision: 1.11 $
%	$Id: kurtosis.m,v 1.11 2003/09/27 09:45:12 schloegl Exp $


if nargin==1,
        DIM=min(find(size(i)>1));
        if isempty(DIM), DIM=1; end;
end;

[R.SUM,R.N,R.SSQ] = sumskipnan(i,DIM);	% sum

R.MEAN 	= R.SUM./R.N;			% mean 
R.SSQ0	= R.SSQ - real(R.SUM).*real(R.MEAN) - imag(R.SUM).*imag(R.MEAN);	% sum square with mean removed

%if flag_implicit_unbiased_estim;    %% ------- unbiased estimates ----------- 
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and SEM are INF
%else
%    n1	= R.N;
%end;

R.VAR  	= R.SSQ0./n1;	     		% variance (unbiased) 
%R.STD	= sqrt(R.VAR);		     	% standard deviation

i       = i - repmat(R.MEAN,size(i)./size(R.MEAN));
%R.CM3 	= sumskipnan(i.^3,DIM)./n1;
R.CM4 	= sumskipnan(i.^4,DIM)./n1;

%R.SKEWNESS = R.CM3./(R.STD.^3);
R = R.CM4./(R.VAR.^2)-3;
