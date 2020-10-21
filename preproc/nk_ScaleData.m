function [x, xmin, xmax, scalefactor] = nk_ScaleData(x, low, high, scalefactor, xmin, xmax, type )
% Usage: [x,xmin,xmax]=scale(x,low,high,type)
% 
% linearly scale x into the range of low to high 
%  it must be that high -low > 0
% if type = 0 (default), low and high are scalars and scaling is uniform in each dim
%         = 1, scaling is done at each individual dimension in the range
%              vectors low and high
% output: x - scaled x,  xmax: original maximum, xmin: original min
% input: x - original x, low: desired min, high: desired max.
% copyright (c) 1996 by Yu hen Hu
% created 9/21/96
% modified 9/22/96: add xmin, xmax to output
% modified 4/5/2001: add type to allow scaling to indivisual ranges in each dim.
% 
if nargin<7,
   type=0;
end

[K,M]=size(x); fl = 0;
range=high-low; % type 0, scalar, type 1, may be a 1 X M vector
if type==0,
   if ~exist('xmin','var') || isempty(xmin), xmin=min(x(:)); end
   if ~exist('xmax','var') || isempty(xmax), xmax=max(x(:)); end
   if ~exist('scalefactor','var') || isempty(scalefactor)
       % Scale to a specific range
       xrange=xmax-xmin;
       scalefactor = range/xrange;
       x = (x - xmin) * scalefactor + low;
   else % Revert scaling to a specific range
       x = (x - low) / scalefactor + xmin;
   end
elseif type==1,
   if size(low,2) == 1, 
       low=ones(1,M)*low; fl = 1;end
   if size(high,2) == 1, 
       high=ones(1,M)*high; fl = 1;end
   if fl, range=high-low; end
   low=diag(diag(low))'; range=diag(diag(range))'; % make them 1 x M vectors
   if ~exist('xmin','var') || isempty(xmin), xmin=min(x); end
   if ~exist('xmax','var') || isempty(xmax), xmax=max(x); end
   if ~exist('scalefactor','var') || isempty(scalefactor)
       xrange=xmax-xmin; % each 1 x M vector
       mask=[xrange<1e-3]; % maskout elements too small, 1 X M
       range=(1-mask).*range+mask; xrange=(1-mask).*xrange+mask; 
       scalefactor = range./xrange;
       x = (x - ones(K,1)*xmin ) * diag(scalefactor) + ones(K,1)*low;
   else
       x = (x - ones(K,1)*low ) / diag(scalefactor) + ones(K,1)*xmin;
   end
end
