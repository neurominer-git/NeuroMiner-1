function [ERR, TBL] = export_scores(handles, batchmode)

warning off

if ~exist('batchmode','var') || isempty(batchmode)
    batchmode = false;
end

switch handles.modeflag
    case 'classification'
        sheetstr = 'Classifier'; 
    case 'regression'
        sheetstr = 'Regressor';  
end

% Do we print the results of the bagged classifier?
prefstr='';
if handles.multi_modal 
    if strcmp(handles.METAstr,'bagged')
        prefstr = '_varB';
    else
        prefstr = sprintf('_var%g', handles.curmodal);
    end
end

% Do we print the results of independent data prediction?
if handles.oocvview
    n=handles.oocvind;
    typestr = sprintf('_A%g_OOCV-%g', handles.curranal, handles.OOCVinfo.Analyses{handles.curranal}.OOCVvec(n));
else
    typestr = sprintf('_A%g_CV', handles.curranal);
end

filename = sprintf('%s%s%s_Predictions', handles.params.TrainParam.SAV.matname, prefstr, typestr);

for i=1:handles.nclass
    sheetname = sprintf('%s%g',sheetstr,i);
    switch handles.modeflag
        case 'classification'
            if ~handles.oocvview
                TBL = handles.BinClass{i}.tbl;
            else
                TBL = handles.OOCV(handles.oocvind).data.tbl(i);
            end
        case 'regression'
            TBL = handles.Regr.tbl;
    end
    [ERR, STATUS, fil, typ] = tbl2file(TBL, filename, sheetname);
    if ~isempty(ERR), break, end
end

% Print multi-group data
if isfield(handles,'MultiClass')
    TBL = handles.MultiClass.tbl;sheetname = 'MultiGroup'; 
    [ERR, STATUS, fil, typ] = tbl2file(TBL, filename, sheetname);
    if handles.oocvview
        if isfield(handles.OOCV(handles.oocvind).data,'tbl_mult')
            sheetname = 'MultiGroup';
            [~, STATUS, fil, typ] = tbl2file(handles.OOCV(handles.oocvind).data.tbl_mult, filename, sheetname);
        end
    end
end

if strcmp(typ,'xls')
   INFO = {'Analysis Index', handles.curranal; ...
            'ID', handles.NM.analysis{handles.curranal}.id; ...
            'Parent Directory', handles.NM.analysis{handles.curranal}.parentdir; ...
            'Root Directory', handles.NM.analysis{handles.curranal}.rootdir; ...
            'ID of NM structure', handles.NM.analysis{handles.curranal}.params.id; ...
            'Time analysis created', handles.NM.analysis{handles.curranal}.meta.TIME; ...
            'User', char(handles.NM.analysis{handles.curranal}.meta.USER); ...
            'OS', sprintf('%s (%s)',handles.NM.analysis{handles.curranal}.meta.OS.name, handles.NM.analysis{handles.curranal}.meta.OS.arch); ...
            'NM version', handles.NM.analysis{handles.curranal}.meta.NM.ver};
    if handles.oocvview
        
        INFO = [INFO;   
            {'OOCV index', handles.oocvind; ...
            'OOCV creation date', handles.OOCVinfo.Analyses{handles.curranal}.date{n}; ...
            'OOCV description', handles.OOCVinfo.Analyses{handles.curranal}.desc{n}; ...
            'OOCV labels provided', handles.OOCVinfo.Analyses{handles.curranal}.labels_known(n); ...
            'OOCV number of subjects', handles.OOCVinfo.Analyses{handles.curranal}.n_subjects_all(n); ...
            'OOCV creation info: System', handles.OOCVinfo.Analyses{handles.curranal}.defs{1}.os_sys; ...
            'OOCV creation info: Version', handles.OOCVinfo.Analyses{handles.curranal}.defs{1}.os_ver; ...
            'OOCV creation info: MATLAB', handles.OOCVinfo.Analyses{handles.curranal}.defs{1}.matlab_ver}];
    end
    xlswrite(fil,INFO,1,'A1');        
end

if STATUS
    if ~batchmode
        msgbox(['Data successfully exported to file ' fil]);
    end
else
    warndlg(['Data NOT successfully exported to file ' fil]);
end

warning on