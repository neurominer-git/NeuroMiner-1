function [ind] = binary2LinearInd(y)

[n,p] = size(y);

for i = 1:n
    ind(i,1) = 0;
    for j = p:-1:1
        if y(i,j) == 1
            ind(i,1) = ind(i,1) + 2^(j-1);
        end
    end
end
ind = ind+1;