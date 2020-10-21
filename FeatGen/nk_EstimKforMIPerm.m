function [MI, alph, optk] = nk_EstimKforMIPerm(Y,Labels,kstart,kend,nperms)

[kX, ix] = size(Y);
kvec = kstart:kend;
kvecX = length(kvec);
miK = zeros(ix,kvecX);
alphK = miK;
meanmiK = zeros(ix,kvecX);
varmiK = zeros(ix,kvecX);
miperm = zeros(nperms,1);

TrainInd = zeros(nperms,kX);
for i=1:nperms
    TrainInd(i,:) = randperm(kX);
end    

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
        miK(i,j) = MIxnyn(Yx,Lb,kvec(j));
        for k=1:nperms
            miperm(k) = MIxnyn(Yx,Lb(TrainInd(i,:)),kvec(j));
        end
        meanmiK(i,j) = mean(miperm); 
        varmiK(i,j) = var(miperm);
        alphK(i,j) = sum(miperm>miK(i,j))/nperms;
    end
end
[minvarK,ind] = min(mean(varmiK));

optk=kvec(ind);
MI=miK(:,ind);
alph=alphK(:,ind);

return