function k=kernelSigmoid(x,y,param)
% sigmoid kernel sigmoid(x,y)=tanh(alpha*(x'*y) + beta)
% Usage:
% k=kernelSigmoid(x,y)
% k=kernelSigmoid(x,y,param)
% x,y: column vectors, or matrices.
% param: scalar, [alpha;beta].
% k, scalar or matrix, the kernel values
% Yifeng Li, May 26, 2011.

k=x'*y;
if nargin<3
   alpha=1;
   beta=0;
else
    alpha=param(1);
    beta=param(2);
end
k=tanh(alpha*k + beta);
end