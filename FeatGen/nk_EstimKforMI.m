function [optk, t] = nk_EstimKforMI(Y,Labels,kstart,kend,K)

optk = [];
[kX, ix] = size(Y);
kvec = kstart:kend;
kvecX = length(kvec);
t = zeros(ix,length(kvec));
mi = zeros(K,1);
miperm = zeros(K,1);

% Generate K-fold subsamples of original data ordering
cv = nk_CVpartition(1,K,Labels);
TrainInd = zeros(length(cv.TrainInd{1}),K);
TrainIndPerm = TrainInd;

% Generate K-fold subsamples of permuted data ordering
PermInd = randperm(kX);
cvperm = nk_CVpartition(1,K,Labels(PermInd));

for i=1:K
    TrainInd(:,i) = cv.TrainInd{i};
    TrainIndPerm(:,i) = cvperm.TrainInd{i};
end
clear cv cvperm

% Loop through all features
for i=1:ix
    Yx = Y(:,i);
    Lb = Labels;
    fprintf('\nFeature %g\t',i)
    if ~sum(any(Y(:,i))), continue, end;
    % Loop through all values of k
    for j=1:kvecX
        fprintf('.')
        % Loop through all CV samples and compute MI
        for k=1:K
            mi(k) = MIxnyn(Yx(TrainInd(:,k)),Lb(TrainInd(:,k)),kvec(j));
            miperm(k) = MIxnyn(Yx(TrainIndPerm(:,k)),Lb(TrainIndPerm(:,k)),kvec(j));
        end
        mmi = mean(mi); vmi = var(mi);
        mmiperm = mean(miperm); vmiperm = var(miperm);
        t(i,j) = (mmi-mmiperm)/sqrt(vmi+vmiperm);
    end
end
[maxt, ind] = max(t);
optk= kvec(ind);

return