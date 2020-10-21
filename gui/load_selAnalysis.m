function handles = load_selAnalysis(handles)

flg1 = 0; flg2=0;  
if isfield(handles,'ExportPredictionsM')
    delete(handles.ExportPredictionsM);
    delete(handles.ExportPerformanceM);
end

for i=1:numel(handles.NM.analysis)
    if handles.NM.analysis{i}.status
        flg1 = flg1+1;
        if isfield(handles.NM.analysis{i},'id')
            popupstr{i} = sprintf('Analysis %g [ %s ]',i, handles.NM.analysis{i}.id);
        else
            popupstr{i} = sprintf('Analysis %g',i);    
        end
        handles.ExportPredictionsM(flg1) = uimenu(handles.ExportPredictionsMenu,'Label',popupstr{i}, 'Tag', num2str(i));
        handles.ExportPerformanceM(flg1) = uimenu(handles.ExportPerformanceMenu,'Label',popupstr{i}, 'Tag', num2str(i));
        if isfield(handles.NM.analysis{i},'visdata')
            flg2=flg2+1;
            if ~isfield(handles,'ExportFeaturesMenu')
                handles.ExportFeaturesMenu = uimenu(handles.ExportMenuItem,'Label','Export Features');
            end
            handles.ExportFeaturesM(flg2) = uimenu(handles.ExportFeaturesMenu,'Label',popupstr{i}, 'Tag', num2str(i));
        end
    else
        popupstr{i} = sprintf('Analysis %g - NOT computed',i);
    end
    handles.n_analyses = flg1;
end

set(handles.selAnalysis, 'String', popupstr);

if flg1>1
     handles.ExportPredictionsM(flg1+1) = uimenu(handles.ExportPredictionsMenu,'Separator', 'on', 'Label','Export predictions of all analyses','Tag','All');
     handles.ExportPerformanceM(flg1+1) = uimenu(handles.ExportPerformanceMenu,'Separator', 'on', 'Label','Export performance metrics of all analyses', 'Tag', 'All');
end

if flg2>1
    handles.ExportFeaturesM(flg2+1) = uimenu(handles.ExportFeaturesMenu, 'Separator', 'on', 'Label','Export feature visualization data of all analyses','Tag','All');
end

if ~isfield(handles,'ExportPredictionsDlg')
    handles.ExportPredictionsDlg = uimenu(handles.ExportMenuItem,'Label','Performance Tabulator', 'Tag', 'PerfTab');
end



