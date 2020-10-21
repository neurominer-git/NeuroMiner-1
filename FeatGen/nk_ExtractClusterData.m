function [meanYC, sdYC] = nk_ExtractClusterData(Y,C)
% function out = nk_ExtractClusterData(Y,C)
% I1: Data matrix with rows = samples, cols = voxels (dimensions)
% I2: Vector indicating cluster membership: length(C) == size(Y,2)!
%
% O1: Output matrix consisting of rows = samples, cols = mean signal of
% sample per cluster
% _________________________________
% (c) Nikolaos Koutsouleris 10/2009

minC = min(C);
maxC = max(C);

nClust = maxC-minC;
nSubj = size(Y,1);

meanYC = zeros(nSubj,nClust);
sdYC = meanYC;

for i=1:nClust
    ind = C==i;
    meanYC(:,i) = mean(Y(:,ind),2);
    %sdYC(:,i) = std(Y(:,C==i),0,2); 
end

return