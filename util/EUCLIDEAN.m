function D = EUCLIDEAN(A, B)

D = sqrt( bsxfun(@plus,sum(A.^2,2),sum(B.^2,2)') - 2*(A*B') );