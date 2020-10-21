function handles = display_main(handles)

h_class         = handles.popupmenu1.Value;
h_classlist     = handles.popupmenu1.String;
h_list          = handles.selModelMeasures.String;
h_val           = handles.selModelMeasures.Value;

switch h_list{h_val}
    
    case 'Classification plot'
     
        handles.pnModelPerf.Visible         = 'off';
        handles.pnVisual.Visible            = 'off';
        handles.pnBinary.Visible            = 'on';
  
        if strcmpi(h_classlist{h_class},'Multi-group classifier')
            handles.oocvview = false;
            handles.cmdExportCobWeb.Visible = 'on';
            handles.selOneVsAll_Info.Enable  = 'on';
            %set(get(handles.pnBinRegrPerfCmd, 'Children'), 'Enable', 'off');
            load_selYAxis(handles)
            load_selModelMeasures(handles)
            
            handles = display_multiclassplot(handles);
            handles = sel_onevsone(handles, handles.selOneVsAll_Info);
            load_selCase(handles,handles.MultiClass.cases)
        else
            %set(get(handles.pnBinRegrPerfCmd, 'Children'), 'Enable', 'on');
            handles.selOneVsAll_Info.Enable = 'off';
            handles.cmdExportCobWeb.Visible = 'off';
            handles.cmdMetricExport.Visible = 'on';
            load_selYAxis(handles)
            load_selModelMeasures(handles)
            
            if strcmp(handles.selCVoocv.Enable,'on') && handles.selCVoocv.Value>1
                handles.oocvview = true;
                handles  = display_classplot_oocv(h_class, handles);
                handles.oocvind = handles.selCVoocv.Value - 1;
                load_selCase(handles,handles.OOCVinfo.Analyses{handles.curranal}.cases{handles.oocvind});
            else
                handles.oocvview = false;
                handles  = display_classplot(h_class, handles);
                load_selCase(handles,handles.BinClass{h_class}.cases)
            end
        end
        
    case 'Regression plot'
        
        handles.pnModelPerf.Visible         = 'off';
        handles.pnVisual.Visible            = 'off';
        handles.pnBinary.Visible            = 'on';
        handles.cmdExportCobWeb.Visible     = 'off';
        handles.cmdMetricExport.Visible     = 'on';
        handles                             = display_regrplot(handles);
        handles                             = binarize_regr(handles);
        load_selCase(handles,handles.Regr.cases)
        
    case 'Visualization results'
      
        handles.pnModelPerf.Visible         ='off';
        handles.pnBinary.Visible            ='off';
        handles.pnVisual.Visible            ='on';
        load_selModality(handles)
        handles = display_visual(handles);
    
    otherwise
      
        handles.pnBinary.Visible            = 'off';
        handles.pnVisual.Visible            = 'off';
        handles.pnModelPerf.Visible         = 'on';
        load_selModelMeasures(handles)
        switch handles.METAstr
            case 'none'
                handles                     = display_modelperf(handles);
            case 'bagged'
                handles                     = display_modelperf_bagged(handles);
        end
end
