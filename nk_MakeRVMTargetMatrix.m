function LabelMatrix = nk_MakeRVMTargetMatrix(LabelVector)

uL = unique(LabelVector);
if numel(uL) > 1
    if min(LabelVector) == -1
        LabelVector(LabelVector == -1) = 2;
        LabelVector(LabelVector == 1) = 1;
    end
    uL(~uL) = [];
    nclass = numel(uL);
else
    nclass = 2;
end

lx = size(LabelVector,1);

LabelMatrix = false(lx,nclass);

for x=1:nclass
    ind = LabelVector == x;
    LabelMatrix(ind,x) = true;
end

return
