function [y] = linearInd2Binary(ind,nLabels)

n = length(ind);

y = -ones(n,nLabels);

for i = 1:n
    bin = dec2bin(ind(i)-1,nLabels);
    for j = 1:nLabels
        if bin(end-j+1) == '1'
           y(i,j) = 1; 
        end
    end
end