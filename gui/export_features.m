function [ERR, TBL] = export_features(handles, batchmode)

warning off

curclass = get(handles.popupmenu1,'Value');
varind = get(handles.selModality,'Value');

if numel(handles.visdata_table(varind).params.NumPred)>1, 
    NumPred = handles.visdata_table(varind).params.NumPred(curclass);
else
    NumPred = handles.visdata_table(varind).params.NumPred;
end

filename = sprintf('%s_A%g_Modality%g_NumPred%g_Cl%g', ...
    handles.params.TrainParam.SAV.matname, ...
    handles.curranal, ...
    handles.visdata_table(varind).params.varind, ...
    NumPred, ...
    curclass);

TBL = handles.visdata_table(varind).tbl(curclass);
ERR = tbl2file(TBL, filename, sprintf('Cl%g',curclass));

if isempty(ERR)
    if ~batchmode
        msgbox(['Feature data successfully exported to file ' filename]);
    end
else
    warndlg(['Feature data NOT successfully exported to file ' filename]);
end

warning on