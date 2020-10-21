function load_selYAxis(handles)

switch handles.modeflag
    case 'classification'
        str = 'classifier scores';
    case 'regression'
        str = 'target predictions';
end

if strcmp(handles.popupmenu1.String{handles.popupmenu1.Value},'Multi-group classifier')
     popuplist = {'Multi-group probabilities derived from CV2-ensembles' , 'Multi-group probabilities derived from similarity averaging'}';
else
    popuplist = {['Mean ' str], ...
                ['Mean ' str ' (95%-CIs)'], ...
                ['Mean ' str ' (SD)']};

    if strcmp(handles.modeflag,'classification')
        popuplist{end+1} = 'Ensemble-based probability scores';
        if isfield(handles.GDdims,'CV2grid')
            popuplist{end+1} = 'Mean ensemble-based probability scores (95%-CIs)';
            popuplist{end+1} = 'Mean ensemble-based probability scores (SD)';
        end
    end
end
set(handles.selYaxis, 'String', popuplist);