function l=log2Inf(n)
% compute log2 of n. log2Inf(0)=1, log2Inf(i)=i if (i~=0).
% n, matrix, vector, or scalar.
% l, the same size as n
% Yifeng Li

n(n==0)=1;
l=log2(n);
end