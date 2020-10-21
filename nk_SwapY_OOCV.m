function NM = nk_SwapY_OOCV(NM)

if ~isfield(NM,'OOCV')
    error('No independent test data found!')
end

fl = questdlg('Are you sure you want to swap the datasets?','Swap datasets','Yes','No','No');
if strcmp(fl,'no'), return; end
    
nV = numel(NM.Y);
nO = numel(NM.OOCV);
if nO>1, 
     selstr = 'select';
     St = 'Independent Test Data Manager'; 
     O = nk_OOCVDataIO('title',St,'list',NM.OOCV,'mode',selstr);
     SelO = O.SelItem;
else
     SelO = 1;
end
OOCV = NM.OOCV{SelO};
nVO = numel(OOCV.Y);
if nVO ~= nV
    error('Independent test data and Training/CV data have unequal number of modalities.')
end
if ~isfield(OOCV,'label') || isempty(OOCV.label)
    error('Independent test data does not contain label information.')
end

NM.OOCV{SelO}.Y = NM.Y;
NM.OOCV{SelO}.cases = NM.cases;
NM.OOCV{SelO}.label = NM.label;
if isfield(NM,'files')
    NM.OOCV{SelO}.files = NM.files;
end
if isfield(NM,'covars')
    NM.OOCV{SelO}.covars = NM.covars;
end
NM.OOCV{SelO}.n_subjects = NM.n_subjects;
NM.OOCV{SelO}.n_subjects_all = NM.n_subjects_all;

NM.Y = OOCV.Y;
NM.cases = OOCV.cases;
NM.label = OOCV.label;
if isfield(OOCV,'files')
    NM.files = OOCV.files;
end
if isfield(NM,'covars')
    NM.covars = OOCV.covars;
end
NM.n_subjects = OOCV.n_subjects;
NM.n_subjects_all = OOCV.n_subjects_all;

NM.OOCV{SelO}.desc = nk_input('Change name of independent test data',0,'s');
NM.OOCV{SelO}.date = datestr(now);
