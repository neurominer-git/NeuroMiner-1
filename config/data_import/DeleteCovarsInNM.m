function NM = DeleteCovarsInNM(NM, covind)
if numel(covind) == size(NM.covars,2);
    NM = rmfield(NM,'covars');
    NM = rmfield(NM,'covnames');
else
    NM.covars(:,covind) = [];
    NM.covnames(covind) = [];
end