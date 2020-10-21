function cv = nk_LOOpartition(Labels)
% =========================================================================
% function cv = nk_LOOpartition(Labels)
% =========================================================================
% (c) Nikolaos Koutsouleris
indn = isfinite(Labels);
K = numel(Labels(indn));
KN = numel(Labels);

trainidxs = cell(1,K); testidxs = cell(1,K); 
ind = 1:KN; cnt=1;

for i=1:KN
    if ~indn(i), continue, end
    indL = true(1,KN); indL(i) = false; 
    trainidxs{cnt} = ind(indL)';
    testidxs{cnt} = ind(~indL)';
    cnt = cnt+1;
end
    
cv.TrainInd = trainidxs;
cv.TestInd = testidxs;

end