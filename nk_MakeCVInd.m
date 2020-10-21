function tInd = nk_MakeCVInd(group, K)

mx = max(group);
N = length(group);
tInd = zeros(N,1);
for g = 1:mx
    h = find(group==g);
    % compute fold id's for every observation in the  group
    q = ceil(K*(1:group(g))/group(g));
    % and permute them to try to balance among all groups
    pq = randperm(K);
    % randomly assign the id's to the observations of this group
    randInd = randperm(group(g));
    tInd(h(randInd))=pq(q);
end