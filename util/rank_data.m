function rnk = rank_data(data,sortmode)

% Sort data
[srt, idxSrt]  = sort(data,sortmode);
% Find where are the repetitions
idxRepeat      = [false diff(srt') == 0];
% Rank with tieds but w/o skipping
rnkNoSkip      = cumsum(~idxRepeat);
% Preallocate rank
rnk            = 1:numel(data);
% Adjust for tieds (and skip)
rnk(idxRepeat) = rnkNoSkip(idxRepeat);
% Sort back
rnk(idxSrt)    = rnk;