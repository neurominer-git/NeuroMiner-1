function distEu=euclidean(A,B)
% compute squared Euclidean distance
% each column is a sample, or observation
[rA,cA]=size(A);
[rB,cB]=size(B);
distEu=nan(cA,cB);
if rA~=rB
    error('Matrix A and B have to have the same number of rows');    
end
for i=1:cA
    distEu(i,:)=sum(( repmat(A(:,i),1,cB) - B).^2);
end
end