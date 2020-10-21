function D = getModelNumDim(h, iy, jy, nP, Pspos, GDFEAT)
    
D=0;
for m = 1 : nP % Parameter combinations
    for k=1:iy % permutations
        for l=1:jy % folds
            D = D + size(GDFEAT{Pspos(m)}{k,l,h},2); 
        end
    end     
end