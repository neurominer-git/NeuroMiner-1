function GDanalysis = nk_MultiEvalPerf(GDanalysis, label, pred, ngroups, fld)

if ~exist('fld','var') || isempty(fld), fld = 'MultiClass'; end
if ~exist('ngroups','var') || isempty(ngroups), ngroups = numel(unique(label(~isnan(label)))); end

ind = ~isnan(pred) & ~isnan(label);
if ~isempty(label)
    errs(ind) = label(ind)~= pred(ind);
    confmatrix = nk_ComputeConfMatrix(label, pred, ngroups);
    % Compute performance measures and assign data to output
    GDanalysis.(fld) = nk_MultiClassAssessConfMatrix(confmatrix, label, pred, errs);
end

