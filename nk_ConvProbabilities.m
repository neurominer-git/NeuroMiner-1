function P = nk_ConvProbabilities(C, ngroups)

m = size(C,1);
P = nan(m,ngroups);
if iscell(C)
    for i=1:m
        if isempty(C{i}) || sum(isnan(C{i}))==size(C{i},2), continue, end
        for j=1:ngroups
            P(i,j) = sum(C{i}==j)/numel(C{i});
        end
    end
else
    indn = ~(sum(isnan(C),2)==size(C,2));
    indnf = ~(sum(isnan(C))==m); 
    for j=1:ngroups
        P(indn,j) = sum(C(indn,indnf)==j,2) / size( C(indn,indnf),2 );
    end
   
end