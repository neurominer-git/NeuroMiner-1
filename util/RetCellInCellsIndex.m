function cols = RetCellInCellsIndex(S, D)

if iscell(D)
    cols = false(1, numel(S));
    for i=1:numel(D)
        ind = strcmp(S,D{i});
        cols(ind) = true;
    end
else
    cols = strcmp(S,D);
end