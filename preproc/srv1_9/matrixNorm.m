function n=matrixNorm(A)
% Calculate the Frobenius norm of matrix A
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 6, 2011

n=sqrt(sum(sum(A.^2)));
end