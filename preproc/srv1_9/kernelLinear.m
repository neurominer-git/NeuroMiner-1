function k=kernelLinear(x,y,param)
% linear kernel
% Usage:
% k=kernelLinear(x,y)
% x,y: column vectors, or matrices.
% param: []
% k, scalar or matrix, the kernel values

if size(x,3)>1
    x=matrizicing(x,3);
    y=matrizicing(y,3);
    x=x';
    y=y';   
end
k=x'*y;
end