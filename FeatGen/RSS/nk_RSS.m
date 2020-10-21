%==========================================================================
%FORMAT feats = nk_RSS(Y, P, transp)
%==========================================================================
%Performs random subspace sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 07 / 2011

function feats = nk_RSS(D, P, transp)

if ~exist('transp','var') || isempty(transp), transp = true; end

nperms          = P.nperms;
nMinFeats       = round((D / 100) * P.nMinFeatsPerc);
nMaxFeats        = round((D / 100) * P.nMaxFeatsPerc);
if nMaxFeats > D, nMaxFeats = D; end
NumFeatsVec     = nMinFeats: nMaxFeats;
NumFeats        = numel(NumFeatsVec);
feats           = zeros(nperms, D);
i=1;
while i <= nperms
    
    % Permute number of features in current subspace
    Pvec = randperm(NumFeats);
    PP = NumFeatsVec(Pvec(1));
    FeatInd = randperm(D);
    FeatVec = FeatInd(1:PP);
    Px = zeros(1, D); Px(FeatVec) = 1;
    
    if any(ismember(Px, feats(1:i,:),'rows'))
        i=i-1; continue
    end
    
    feats(i, : ) = Px;
    i=i+1;
end

if transp, feats = feats'; end