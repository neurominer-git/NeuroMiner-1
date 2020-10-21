function k=kernelPoly(x,y,param,varargin)
% Polynomial kernel k=(Gamma.*(x'*y)+ Coefficient).^Degree;
% Usage:
% k=kernelPoly(x,y)
% k=kernelPoly(x,y,param)
% x,y: column vectors, or matrices.
% param: [Gamma;Coefficient;Degree]
% k, scalar or matrix, the kernel values
% Yifeng Li, May 26, 2011.

if nargin<3
    Gamma=1;
    Coefficient=0;
    Degree=2;
else
    Gamma=param(1);
    Coefficient=param(2);
    Degree=param(3);
end
k=(Gamma.*(x'*y)+ Coefficient).^Degree;
end