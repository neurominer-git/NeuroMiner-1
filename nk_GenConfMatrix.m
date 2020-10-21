function [confmatrix, err] = nk_GenConfMatrix(labels, pred)

ml = max(labels); lx = length(labels); confmatrix = zeros(ml,ml);
err = labels ~= pred;

for i=1:lx
    confmatrix(labels(i),pred(i)) = confmatrix(labels(i),pred(i)) + 1;
end

return