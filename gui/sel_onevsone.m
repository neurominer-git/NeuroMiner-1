function handles = sel_onevsone(handles, hObject)

rowind = get(hObject,'Value');
switch rowind
    case 1
        %% Display ROC
        delete(findall(handles.figure1,'Tag','AnnotPerfMeas'))
        [handles.hroc, handles.hroc_random] = display_roc(handles);
        %% Display Cobweb
        [handles.hspider, handles.MultiClass.misclass_confusion] = nk_PlotCobWeb(handles.MultiClass.confusion_matrix, handles.NM.groupnames, handles.axes5);
        %handles.axes4.Visible='off'; handles.axes3.Visible='off';
        %handles.axes4.Title.Visible = 'off'; handles.axes4.Title.Visible = 'off';
        cla(handles.axes4);title(handles.axes4,{''});
        cla(handles.axes3);title(handles.axes3,{''});
        handles.selOneVsAll_Info.Visible = 'on';
        handles.txtPretestProb.Visible = 'off';
        handles.cmdExportPies.Visible = 'off';
        handles.cmdExportCobWeb.Visible = 'on';
        handles.cmdMetricExport.Visible = 'off';
        
    otherwise
        if isfield(handles,'hspider'),handles.hspider.Title.Visible='off'; end
        handles.cmdExportPies.Visible = 'on';
        handles.cmdMetricExport.Visible = 'on';
        handles.cmdExportCobWeb.Visible = 'off';
        
        %% Display ROC
        %handles.txtPerf.Visible='on';
        legend off %set(handles.txtPerf,'visible','on');
        
        [handles.hroc, handles.hroc_random] = display_roc(handles);
        
        %% Display contingency info
        handles.h_contiginfo = display_contigmat(handles);
        
        %% Display pie charts
        [handles.h1pie, handles.h2pie] = display_piecharts(handles);
        
end
%% Display confusion matrix
handles.h_contig = display_contigplot(handles);