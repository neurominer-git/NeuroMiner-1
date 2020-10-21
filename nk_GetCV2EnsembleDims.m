function Ydims = nk_GetCV2EnsembleDims(D)

[nperms, nfolds, nclass] = size(D);

Ydims = zeros(nclass,1);
%Xdims = ones(nclass,1);

for i=1:nperms
    for j=1:nfolds
        for curclass=1:nclass
            DY = size(D{i,j,curclass},2);
            Ydims(curclass) = Ydims(curclass) + DY;
        end
    end
end
