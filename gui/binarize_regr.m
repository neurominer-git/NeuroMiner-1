function handles = binarize_regr(handles)

label = handles.NM.label(:,handles.curlabel); 
indnan = ~isnan(label); label = label(indnan); 
pred = handles.Regr.mean_predictions(indnan);
m = str2double(get(handles.txtBinarize,'String')); 
if isempty(m)
    errordlg('Enter numeric threshold')    
elseif m <= min(label) || m >= max(label)
    errordlg(sprintf('Threshold out of target range [%g %g]',min(label),max(label)));
else
    handles.Regr.b_label = label; handles.Regr.b_label(label>=m) = 1; handles.Regr.b_label(label<m) = -1;
    handles.Regr.b_pred = pred; 
    handles.Regr.b_pred = handles.Regr.b_pred - m;
    [handles.Regr.X, ...
     handles.Regr.Y, ...
     handles.Regr.T, ...
     handles.Regr.AUC] = perfcurve2(handles.Regr.b_label, handles.Regr.b_pred, 1);
    handles.Regr.contigmat = ALLPARAM(handles.Regr.b_label, handles.Regr.b_pred);
    handles.Regr.contigmat.BINARIZATION_THRESHOLD = m;
    handles = display_contigmat(handles);
end

%% Display ROC
[handles.hroc, handles.hroc_random] = display_roc(handles);

%% Display pie charts
[handles.h1pie, handles.h2pie] = display_piecharts(handles);