function C=nm_nanplus(A,B);
% Adding Two Matrices, but ingnoreing Nan-entries
A(isnan(A))=0;
B(isnan(B))=0;
C=A+B;