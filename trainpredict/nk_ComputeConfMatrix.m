function confmatrix = nk_ComputeConfMatrix(label, pred, ngroups)

ind = ~isnan(pred) & ~isnan(label);
confmatrix = zeros(ngroups);
lx = size(label,1);
% Compute confusion matrix
for i=1:lx
    if ~ind(i), continue, end
    confmatrix(label(i),pred(i)) = confmatrix(label(i),pred(i)) + 1;
end
   
