function R = nk_FScore(Y, L)

indP = L == 1; 
indM = L == -1;
nYp = sum(indP); 
nYm = sum(indM);
mY = mean(Y); 
mYp = mean(Y(indP,:)); 
mYm = mean(Y(indM,:));

R = ((mYp - mY).^2 + (mYm - mY).^2) / ( sum((Y(indP,:)-mYp) ^ 2) / (nYp - 1) + sum((Y(indM,:)-mYm) ^ 2) / (nYm - 1) ); 

end