function C = AUtoPhysicalUnits(X,B)

% INPUT
% X: data matrix (DxN; D: dimensionality of input data, N: number of samples)
% B: estimated components (DxK; K: number of components)
%
% Note that this function can also be used to calculate coefficients for out of sample subjects 
%
% OUTPUT
% C: new matrix (KxN) - values in physical units

    
% normalize to sum to 1
Blen = sum(B,1);
if any(Blen==0)
	Blen(Blen==0) = 1;
end
nB = bsxfun(@times,B,1./Blen) ;
    
C = nB' * X ;