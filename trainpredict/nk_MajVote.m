function [probmax, prob, randvote] = nk_MajVote(H, M, curlabel)

if ~exist('curlabel','var'), curlabel=1; end
uM      = unique(M,'stable'); 
nuM     = numel(uM);
l       = size(H,1);
if iscell(H), H = cell2mat(H(:,curlabel)); end
if isequal(M, [1 -1]), H = sign(H);end
prob = zeros(l,nuM); probmax=zeros(l,1); randvote = false(l,1);
for i=1:l
    [prob(i,:), probmax(i) , randvote(i)] = MajVote(H(i,:), uM, nuM);
end
if sum(randvote)
    warning('NeuroMiner:MajVote: %g tied vote counts encountered!', sum(randvote))
end
% Check whether you have to through coin

function [prob, votes, randvote] = MajVote(H,uM, nuM)

randvote = false;
l       = nm_nancount(H);
prob    = arrayfun( @(i) rdivide(nm_nansum(H==uM(i)),l),1:nuM);
[probmax,ind] = max(prob);
% Check whether a coin has to be thrown
lu = find(prob==probmax); if numel(lu)>1, ind = randperm(numel(lu),1); ind = lu(ind); randvote = true; end
votes = uM(ind);