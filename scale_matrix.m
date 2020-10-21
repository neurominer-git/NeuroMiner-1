function [vscaled,minv,maxv]= scale_matrix(v,minv,maxv)
% scale_matrix scales a matrix into [0 1] along the vertical direction,
% this is a simplified version
%
% SYNTAX
%
% Inputs: v: m*n matrix;
% Outputs:
% Wenbin, 08-Mar-2007

vr =size(v,1);
if vr <2
error('The rows of input matrix must be greater than 1');
end
if nargin < 2
    maxv=max(v);
    minv=min(v);
end

diff =maxv-minv;
diff(diff ==0) =1;% Avoid dividing by zero
vscaled =(v-repmat(minv, vr,1))./(repmat(diff,vr,1));