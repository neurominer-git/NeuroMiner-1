function [val, perms] = rfe_multi_perm_featcombs(D, T, L, nclass, nperms, ngroups)
% =========================================================================
% FORMAT [val, perms] = rfe_multi_perm_featcombs(D, T, L, nclass, nperms)
% =========================================================================
% Identify fittest combination of features across dichotomizers using a
% permutation approach
% =========================================================================
% Nikolaos Koutsouleris, 01/2019

nfeats = size(D,3);
perms = zeros(nclass, nfeats, nperms);
for curclass=1:nclass
    for curperm=1:nperms
        perms(curclass,:,curperm) = randperm(nfeats);
    end
end
val = zeros(nfeats, nperms);

for i=1:nperms
    tD = D; tT = T;
    for curclass=1:nclass
        tD(:,curclass,:) = D(:,curclass, perms(curclass,:,i));
        tT(:,curclass,:) = T(:,curclass, perms(curclass,:,i));
    end
    for j=1:nfeats
        val(j,i) = nk_MultiEnsPerf(tD(:,:,j), tT(:,:,j), L, 1:nclass, ngroups);
    end
end
