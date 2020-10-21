function MutualInformation = nk_MutualInfo(Y, labels)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% function MutualInformation = nk_MutualInfo(symY, labels)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ix = size(Y,2);
MutualInformation = zeros(ix,1);
try
    parfor i=1:ix
        sY = Y(:,i);
        MutualInformation(i) = mutualinfo(sY,labels);
    end
catch
    for i=1:ix
        MutualInformation(i) = mutualinfo(Y(:,i),labels);
    end
end
return