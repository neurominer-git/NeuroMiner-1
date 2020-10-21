function X=l1NNLSKernelBatch(A,B,option)
% Newton's algorithm to solve the kernel l1LS problem:
% min f(X)=1/2||phi(B)-phi(A)X||_2^2 + lambda||X||_1 where X is a sparse matrix
% Yifeng Li
% Feb. 08, 2012

if nargin<3
   option=[]; 
end
optionDefault.lambda=0.1;
optionDefault.kernel='rbf';
optionDefault.param=2^0;
option=mergeOption(option,optionDefault);

% compute kernel matrices
n=size(A,2);
numB=size(B,2);
AtA=computeKernelMatrix(A,A,option);
AtB=computeKernelMatrix(A,B,option);
BtB=zeros(numB,1);
for i=1:numB
    BtB(i)=computeKernelMatrix(B(:,i),B(:,i),option);
end
% % normalize kernel
if strcmp(option.kernel,'randrbf')
    [AtA,~,BtB,~,AtB]=normalizeKernelMatrix(AtA,diag(BtB),AtB);
    BtB=diag(BtB);
end
% compute coefficients
X=nan(n,numB);
AtAInv=eye(n)/(AtA+2^-64*eye(n));
for i=1:numB
    X(:,i)=l1NNLSKernel(AtA,AtB(:,i),BtB(i),option.lambda,AtAInv,option);
end
end