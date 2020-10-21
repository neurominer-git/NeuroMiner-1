function [ERR, TBL] = export_performance(handles, batchmode)

warning off 

switch handles.modeflag
    case 'classification'
        sheetstr = 'Classifier'; 
    case 'regression'
        sheetstr = 'Regressor';  
end

filename = sprintf('%s_A%g_PredictionMetrics', ...
        handles.params.TrainParam.SAV.matname, handles.curranal);

for i=1:handles.nclass
    sheetname = sprintf('%s%g',sheetstr,i);
    switch handles.modeflag
        case 'classification'
            TBL = handles.BinClass{i}.tbl_cont;
        case 'regression'
            TBL = handles.Regr.tbl_cont;
            if isfield(handles.Regr,'contigmat')
                tbl_cont.rownames   = fieldnames(handles.Regr.contigmat);
                tbl_cont.array      = struct2cell( handles.Regr.contigmat);
                remind = find(strcmp('FPRvec',tbl_cont.rownames) | strcmp('TPRvec', tbl_cont.rownames) | strcmp('X',tbl_cont.rownames));
                tbl_cont.array(remind) = [];
                tbl_cont.array = cell2mat(tbl_cont.array);
                tbl_cont.rownames(remind) = [];
                TBL.array = [ TBL.array; tbl_cont.array ];
                TBL.rownames = [TBL.rownames; tbl_cont.rownames ];
            end
    end
    ERR = tbl2file(TBL, filename, sheetname);
end
if isfield(handles,'MultiClass')
    sheetname = sprintf('MultiClass');
    ERR = tbl2file(handles.MultiClass.tbl_cont, filename, sheetname);
end

warning on