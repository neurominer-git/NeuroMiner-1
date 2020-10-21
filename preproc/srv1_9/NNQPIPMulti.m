function X=NNQPIPMulti(AtA,AtB,BtB,option)
% Newton's algorithm to solve the kernel l1LS problem:
% min f(X)=1/2||phi(B)-phi(A)X||_2^2 + lambda||X||_1 where X is a sparse matrix
% Yifeng Li
% Feb. 08, 2012

if nargin<3
   option=[]; 
end
optionDefault=[];
option=mergeOption(option,optionDefault);

% compute kernel matrices
n=size(AtA,1);
numB=size(AtB,2);

% compute coefficients
X=nan(n,numB);
AtAInv=eye(n)/(AtA+2^-64*eye(n));
for i=1:numB
    X(:,i)=NNQPIP(AtA,AtB(:,i),BtB(i),AtAInv,option);
end
end