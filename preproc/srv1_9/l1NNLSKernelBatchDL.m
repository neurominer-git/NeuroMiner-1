function X=l1NNLSKernelBatchDL(AtA,AtB,BtB,option)
% Newton's algorithm to solve the kernel l1LS problem:
% min f(X)=1/2||phi(B)-phi(A)X||_2^2 + lambda||X||_1 where X is a sparse matrix
% Yifeng Li
% Feb. 08, 2012

if nargin<3
   option=[]; 
end
optionDefault.lambda=0.1;
option=mergeOption(option,optionDefault);

% compute kernel matrices
n=size(AtA,1);
numB=size(AtB,2);

% compute coefficients
X=nan(n,numB);
AtAInv=eye(n)/(AtA+2^-64*eye(n));
for i=1:numB
    X(:,i)=l1NNLSKernel(AtA,AtB(:,i),BtB(i),option.lambda,AtAInv,option);
end
end