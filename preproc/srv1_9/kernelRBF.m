function k=kernelRBF(x,y,param)
% rbf kernel rbf(x,y)=exp((-(1/sigma^2).*(|x-y|.^2));
% Usage:
% k=kernelRBF(x,y)
% k=kernelRBF(x,y,param)
% x,y: column vectors, or matrices.
% param: scalar, [sigma].
% k, scalar or matrix, the kernel values
% Yifeng Li, May 26, 2011.


if nargin<3
    sigma=1;
else
    sigma=param(1,1);
    if sigma==0
        error('sigma must not be zero!');
    end
end
k=exp((-(1/sigma^2)).*(repmat(sum(x.^2,1)',1,size(y,2))-2*(x'*y)+repmat(sum(y.^2,1),size(x,2),1)));
end