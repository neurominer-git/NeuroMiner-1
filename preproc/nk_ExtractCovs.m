function [ecovs, eind] = nk_ExtractCovs(covs, ACTPARAM)

lact = numel(ACTPARAM); ecovs=[];eind=[];

for i=1:lact
    if ~isempty(ACTPARAM{i}) && strcmp(ACTPARAM{i}.cmd,'correctnuis')
        ecovs = covs(:,ACTPARAM{i}.COVAR); break
    end
   
end

for i=1:lact
    if ~isempty(ACTPARAM{i}) && strcmp(ACTPARAM{i}.cmd,'normalize')
        eind = covs(:,ACTPARAM{i}.IND); break
    end
end