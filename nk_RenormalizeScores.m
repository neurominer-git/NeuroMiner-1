function xDs = nk_RenormalizeScores(Dr,Ds)

nR = size(Dr,1); nS = size(Ds,1);
X = floor(nR/nS);

[sDr, indR] = sort(Dr,'descend');
[~, indS] = sort(Ds,'descend');

xDs = zeros(nS,1);

vec = 1: X : nR;
for i = 1:numel(vec)
   
    xDs(i) = mean(sDr(vec(i):vec(i+1))); 
    
end