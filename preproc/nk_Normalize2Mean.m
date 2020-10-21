function sY = nk_Normalize2Mean(Y)

[n, m] = size(Y);
sY = zeros(n,m);

for i=1:n
    
    fprintf('\nNormalizing subject %g',i)
    
    Yi = Y(i,:);
    indNonZero = ~(Yi == 0);
    meani = median(Yi(indNonZero));
    sY(i,indNonZero) = Yi(indNonZero) * 100 / meani;  
    
end