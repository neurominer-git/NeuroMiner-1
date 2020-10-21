function NM = nk_Move2OOCV(NM,I)

if isfield(NM,'OOCV'),
    OOCVind  = numel(NM.OOCV)+1;
else
    OOCVind = 1;
end

NM.OOCV{OOCVind}.desc   = 'OOCV created by extracting cases from CV';
NM.OOCV{OOCVind}.date   = date;
NM.OOCV{OOCVind}.ind    = I;
NM.OOCV{OOCVind}.fldnam = 'OOCV';

nM = numel(NM.Y);

NM.OOCV{OOCVind}.cases = NM.cases(I);
NM.OOCV{OOCVind}.label = NM.label(I,:);
if isfield(NM,'covars')
    NM.OOCV{OOCVind}.covars = NM.covars(I,:);
    NM.OOCV{OOCVind}.covnames = NM.covnames;
end

NM.OOCV{OOCVind}.groupnames = NM.groupnames;
if strcmp(NM.modeflag,'classification')
    uL = unique(NM.label);
    uL(isnan(uL))=[];
    for i=1:numel(uL)
        NM.OOCV{OOCVind}.n_subjects(i) = sum(NM.label(I)==uL(i));
    end
else
    NM.OOCV{OOCVind}.n_subjects = sum(I);
end
NM.OOCV{OOCVind}.n_subjects_all = sum(I);

for i=1:nM
    NM.OOCV{OOCVind}.Y{i} = NM.Y{i}(I,:);
    NM.Y{i}(I,:) = [];
    if ~isempty(NM.files{i})
        NM.OOCV{OOCVind}.files{i} = NM.files{i}(I,:);
        NM.files{i}(I,:)=[];
    end
end

NM.cases(I)=[];
NM.covars(I,:)=[];
NM.label(I,:)=[];

if strcmp(NM.modeflag,'classification')
    uL = unique(NM.label);
    uL(isnan(uL))=[];
    for i=1:numel(uL)
        NM.n_subjects(i) = sum(NM.label==uL(i));
    end
else
    NM.n_subjects = sum(~I);
end

NM.n_subjects_all = sum(~I);

if ~isempty(NM.TrainParam.RAND.CV2LCO)
    NM.TrainParam.RAND.CV2LCO.ind(I)=[];
end


