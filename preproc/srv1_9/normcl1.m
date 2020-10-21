function Y=normcl1(Y)
% normalization, each column of Y has l1 norm 1
% Y is nonnegative matrix

Y=abs(Y);
sumY=sum(Y,1);
[r,c]=size(Y);
Y=Y./(repmat(sumY,r,1));
end
